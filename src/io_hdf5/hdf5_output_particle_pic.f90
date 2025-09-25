!=================================================================================================================================
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

MODULE MOD_HDF5_Output_Particles_PIC
#if defined(PARTICLES)
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
#if !(USE_FV) || (USE_HDG)
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_IO_HDF5
USE MOD_HDF5_output
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: WriteNodeSourceExtToHDF5
PUBLIC :: WriteElectroMagneticPICFieldToHDF5
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteNodeSourceExtToHDF5(OutputTime)
!===================================================================================================================================
! Write NodeSourceExt (external charge density) field to HDF5 file
!===================================================================================================================================
! MODULES
USE MOD_io_HDF5
USE MOD_Globals
USE MOD_PreProc
USE MOD_Dielectric_Vars    ,ONLY: NodeSourceExtGlobal
USE MOD_Mesh_Vars          ,ONLY: MeshFile,offsetElem,nElems
USE MOD_Mesh_Tools         ,ONLY: GetCNElemID
USE MOD_Globals_Vars       ,ONLY: ProjectName
USE MOD_PICDepo_Vars       ,ONLY: NodeSourceExt,NodeVolume,DoDeposition
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_Particle_Mesh_Vars ,ONLY: ElemNodeID_Shared,NodeInfo_Shared,nUniqueGlobalNodes
USE MOD_TimeDisc_Vars      ,ONLY: iter
USE MOD_Interpolation_Vars ,ONLY: NodeType,NodeTypeVISU,Nmin,Nmax
USE MOD_Interpolation      ,ONLY: GetVandermonde
USE MOD_DG_vars            ,ONLY: N_DG_Mapping,nDofsMapping
#if USE_MPI
USE MOD_PICDepo            ,ONLY: ExchangeNodeSourceExtTmp
#endif /*USE_MPI*/
USE MOD_HDF5_Output_ElemData,ONLY: WriteAdditionalElemData
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: OutputTime
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,PARAMETER              :: nVarOut=1
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
CHARACTER(LEN=255)             :: FileName,DataSetName
INTEGER                        :: iElem,iMax,CNElemID
REAL                           :: NodeSourceExtEqui(1:nVarOut,0:1,0:1,0:1),sNodeVol(1:8)
INTEGER                        :: NodeID(1:8)
! p-adaption output
TYPE VdmType
  REAL,ALLOCATABLE :: Vdm(:,:)                              !< Vandermonde mapping from equidistant (visu) to NodeType node set
END TYPE VdmType

TYPE(VdmType), DIMENSION(:), ALLOCATABLE :: Vdm_EQ_N        !< Array to store all Vandermonde matrices depending on Nloc

TYPE NSEType
  REAL,ALLOCATABLE :: U(:,:,:,:)                            !< NodeSourceExtEqui mapped to Nloc
END TYPE NSEType

TYPE(NSEType), DIMENSION(:), ALLOCATABLE :: NodeSourceExt_N !< Array to store all NodeSourceExtEqui depending on Nloc

REAL,ALLOCATABLE               :: U_N_2D_local(:,:)
INTEGER                        :: iDOF, nDOFOutput, offsetDOF, Nloc, i, j, k
!===================================================================================================================================
ALLOCATE(StrVarNames(1:nVarOut))
StrVarNames(1)='NodeSourceExt'

! build necessary mappings
ALLOCATE(Vdm_EQ_N(Nmin:Nmax))
ALLOCATE(NodeSourceExt_N(Nmin:Nmax))
DO Nloc = Nmin, Nmax
  ALLOCATE(Vdm_EQ_N(Nloc)%Vdm(0:Nloc,0:1))
  CALL GetVandermonde(1, NodeTypeVISU, Nloc, NodeType, Vdm_EQ_N(Nloc)%Vdm(0:Nloc,0:1), modal=.FALSE.)
  ALLOCATE(NodeSourceExt_N(Nloc)%U(1:1,0:Nloc,0:Nloc,0:Nloc))
END DO ! Nloc = Nmin, Nmax

! Preparing U_N_2D_local array for output as DG_Solution
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

! Skip MPI communication in the first step as nothing has been deposited yet
IF(iter.NE.0)THEN

#if USE_MPI
! Communicate the NodeSourceExtTmp values of the last boundary interaction before the state is written to .h5
! Only call when deposition is active (otherwise this routine only writes the old array from the restart file to keep the data)
IF(DoDeposition) CALL ExchangeNodeSourceExtTmp()
#endif /*USE_MPI*/

end if ! iter.NE.0

! Loop over all elements and store charge density values in equidistantly distributed nodes of PP_N=1
sNodeVol(1:8) = 1. ! default
! This array is not allocated when DoDeposition=F, however, the previously calculated surface charge might still be required in
! the future, when DoDeposition is activated again. Therefore, read the old data and store in the new state file.
IF(.NOT.ALLOCATED(NodeSourceExt)) THEN
  ALLOCATE(NodeSourceExt(1:nUniqueGlobalNodes))
  NodeSourceExt = 0.
END IF
! Write into 2D array
iDOF = 0
DO iElem=1,nElems
  ! Copy values to equidistant distribution
  CNElemID = GetCNElemID(iElem+offsetElem)
  NodeID = NodeInfo_Shared(ElemNodeID_Shared(:,CNElemID))
  IF(DoDeposition) sNodeVol(1:8) = 1./NodeVolume(NodeID(1:8)) ! only required when actual deposition is performed
  NodeSourceExtEqui(1,0,0,0) = NodeSourceExt(NodeID(1))*sNodeVol(1)
  NodeSourceExtEqui(1,1,0,0) = NodeSourceExt(NodeID(2))*sNodeVol(2)
  NodeSourceExtEqui(1,1,1,0) = NodeSourceExt(NodeID(3))*sNodeVol(3)
  NodeSourceExtEqui(1,0,1,0) = NodeSourceExt(NodeID(4))*sNodeVol(4)
  NodeSourceExtEqui(1,0,0,1) = NodeSourceExt(NodeID(5))*sNodeVol(5)
  NodeSourceExtEqui(1,1,0,1) = NodeSourceExt(NodeID(6))*sNodeVol(6)
  NodeSourceExtEqui(1,1,1,1) = NodeSourceExt(NodeID(7))*sNodeVol(7)
  NodeSourceExtEqui(1,0,1,1) = NodeSourceExt(NodeID(8))*sNodeVol(8)

  ! Map equidistant distribution to G/GL (current node type)
  Nloc = N_DG_Mapping(2,iElem+offsetElem)
  CALL ChangeBasis3D(1, 1, Nloc, Vdm_EQ_N(Nloc)%Vdm, NodeSourceExtEqui(1:1,0:1   ,0:1   ,0:1)   , &
                                               NodeSourceExt_N(Nloc)%U(1:1,0:Nloc,0:Nloc,0:Nloc))
  DO k=0,Nloc; DO j=0,Nloc; DO i=0,Nloc
    iDOF = iDOF + 1
    U_N_2D_local(1:nVarOut,iDOF)   = NodeSourceExt_N(Nloc)%U(1:1,i,j,k)
  END DO; END DO; END DO
END DO!iElem

DEALLOCATE(Vdm_EQ_N)
DEALLOCATE(NodeSourceExt_N)

! Write data twice to .h5 file
! 1. to _State_.h5 file (or restart)
! 2. to separate file (for visu)
#if USE_DEBUG
iMax=2 ! write to state and to a separate file (for debugging)
#else
iMax=1 ! write to state file
#endif /*USE_DEBUG*/
DO i = 1, iMax
  IF(i.EQ.1)THEN
    ! Write field to _State_.h5 file (or restart)
    FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_State',OutputTime))//'.h5'
    DataSetName='DG_SourceExt'
  ELSE
    ! Generate skeleton for the file with all relevant data on a single processor (MPIRoot)
    ! Write field to separate file for debugging purposes
    CALL GenerateFileSkeleton('NodeSourceExtGlobal',nVarOut,StrVarNames,TRIM(MeshFile),OutputTime,FileNameOut=FileName)
#if USE_MPI
    CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif
    IF(MPIRoot)THEN
      CALL OpenDataFile(FileName,create=.FALSE.,single=.TRUE.,readOnly=.FALSE.,communicatorOpt=MPI_COMM_PICLAS)
      CALL WriteAttributeToHDF5(File_ID,'VarNamesNodeSourceExtGlobal',nVarOut,StrArray=StrVarNames)
      CALL CloseDataFile()
    END IF ! MPIRoot
    DataSetName='DG_Solution'

    ! Write 'Nloc' array to the .h5 file, which is required for 2D DG_Solution conversion in piclas2vtk
    CALL WriteAdditionalElemData(FileName,ElementOutNloc)
  END IF ! i.EQ.2

  ! Associate construct for integer KIND=8 possibility
  ASSOCIATE(nVarOut         => INT(nVarOut,IK)           ,&
            nDofsMapping    => INT(nDofsMapping,IK)      ,&
            nDOFOutput      => INT(nDOFOutput,IK)        ,&
            offsetDOF       => INT(offsetDOF,IK)         )
    CALL GatheredWriteArray(FileName,create=.FALSE.,&
                          DataSetName = TRIM(DataSetName) , rank = 2                , &
                          nValGlobal  = (/nVarOut         , nDofsMapping/)          , &
                          nVal        = (/nVarOut         , nDOFOutput/)            , &
                          offset      = (/0_IK            , offsetDOF/)             , &
                          collective  = .TRUE.            , RealArray = U_N_2D_local)
  END ASSOCIATE
  SDEALLOCATE(U_N_2D_local)
END DO ! i = 1, 2

SDEALLOCATE(NodeSourceExtGlobal)
SDEALLOCATE(StrVarNames)
END SUBROUTINE WriteNodeSourceExtToHDF5


!===================================================================================================================================
!> Store the magnetic filed acting on particles at each DOF for all elements to .h5
!===================================================================================================================================
SUBROUTINE WriteElectroMagneticPICFieldToHDF5()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Globals_Vars           ,ONLY: ProjectName
USE MOD_Mesh_Vars              ,ONLY: offsetElem,nElems,MeshFile,N_VolMesh
USE MOD_Output_Vars            ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Interpolation_Vars     ,ONLY: NodeType
USE MOD_PICInterpolation_tools ,ONLY: GetExternalFieldAtParticle,GetEMField
USE MOD_Restart_Vars           ,ONLY: RestartTime
USE MOD_Interpolation_Vars     ,ONLY: N_Inter,Nmax,Nmin
USE MOD_DG_Vars                ,ONLY: nDofsMapping,N_DG_Mapping
USE MOD_HDF5_Output_ElemData   ,ONLY: WriteAdditionalElemData
USE MOD_IO_HDF5                ,ONLY: ElementOutNloc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: FileName
CHARACTER(LEN=255),ALLOCATABLE :: StrVarNames(:)
INTEGER,PARAMETER              :: nVarOut=6
REAL,ALLOCATABLE               :: U_N_2D_local(:,:)
REAL                           :: StartT,EndT
INTEGER                        :: iElem,k,i,j,iDOF,nDOFOutput,offsetDOF,Nloc
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO')' WRITE PIC EM-FIELD TO HDF5 FILE...'
GETTIME(StartT)

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

DO iElem=1,PP_nElems
  Nloc = N_DG_Mapping(2,iElem+offsetElem)
  DO k=0,Nloc
    DO j=0,Nloc
      DO i=0,Nloc
        iDOF = iDOF + 1
        ASSOCIATE( x => N_VolMesh(iElem)%Elem_xGP(1,i,j,k), y => N_VolMesh(iElem)%Elem_xGP(2,i,j,k), z => N_VolMesh(iElem)%Elem_xGP(3,i,j,k))
          ! Superposition of the external and calculated electromagnetic field
          !   GetExternalFieldAtParticle : Get the 1 of 4 external fields (analytic, variable, etc.) at position x,y,z
          !                   GetEMField : Evaluate the electro-(magnetic) field using the reference position and return the field
          U_N_2D_local(1:nVarOut,iDOF) = GetExternalFieldAtParticle((/x,y,z/)) &
                                       + GetEMField(iElem,(/N_Inter(Nloc)%xGP(i),N_Inter(Nloc)%xGP(j),N_Inter(Nloc)%xGP(k)/))
        END ASSOCIATE
      END DO ! i
    END DO ! j
  END DO ! k
END DO ! iElem=1,PP_nElems


! Create dataset attribute "VarNames"
ALLOCATE(StrVarNames(1:nVarOut))
StrVarNames(1)='ElectricFieldX'
StrVarNames(2)='ElectricFieldY'
StrVarNames(3)='ElectricFieldZ'
StrVarNames(4)='MagneticFieldX'
StrVarNames(5)='MagneticFieldY'
StrVarNames(6)='MagneticFieldZ'

! Generate skeleton for the file with all relevant data on a single proc (MPIRoot)
FileName=TRIM(ProjectName)//'_PIC-EMField.h5'
IF(MPIRoot) THEN
  CALL OpenDataFile(TRIM(FileName),create=.TRUE.,single=.TRUE.,readOnly=.FALSE.,userblockSize=userblock_total_len)
  ! Write file header
  CALL WriteHDF5Header('BField',File_ID) ! File_Type='BField'
  ! Write dataset properties "Time","MeshFile","NextFile","NodeType","VarNames"
  CALL WriteAttributeToHDF5(File_ID,'N',1,IntegerScalar=PP_N)
  CALL WriteAttributeToHDF5(File_ID , 'Nmin'     , 1       , IntegerScalar = Nmin)
  CALL WriteAttributeToHDF5(File_ID , 'Nmax'     , 1       , IntegerScalar = Nmax)
  CALL WriteAttributeToHDF5(File_ID,'MeshFile',1,StrScalar=(/TRIM(MeshFile)/))
  CALL WriteAttributeToHDF5(File_ID,'NodeType',1,StrScalar=(/NodeType/))
  CALL WriteAttributeToHDF5(File_ID , 'VarNames' , nVarOut , StrArray      = StrVarNames)
  CALL WriteAttributeToHDF5(File_ID,'Time'    ,1,RealScalar=RestartTime)
  CALL CloseDataFile()
  ! Add userblock to hdf5-file
  CALL copy_userblock(TRIM(FileName)//C_NULL_CHAR,TRIM(UserblockTmpFile)//C_NULL_CHAR)
END IF
#if USE_MPI
CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
#endif /*USE_MPI*/

! Write 'Nloc' array to the .h5 file, which is required for 2D DG_Solution conversion in piclas2vtk
CALL WriteAdditionalElemData(FileName,ElementOutNloc)

! Associate construct for integer KIND=8 possibility
ASSOCIATE(nVarOut         => INT(nVarOut,IK)      ,&
          nDofsMapping    => INT(nDofsMapping,IK) ,&
          nDOFOutput      => INT(nDOFOutput,IK)   ,&
          offsetDOF       => INT(offsetDOF,IK)    )
CALL GatheredWriteArray(FileName,create = .FALSE.                  , &
                           DataSetName  = 'DG_Solution', rank = 2  , &
                           nValGlobal   = (/nVarOut, nDofsMapping/), &
                           nVal         = (/nVarOut, nDOFOutput/)  , &
                           offset       = (/0_IK   , offsetDOF/)   , &
                           collective   = .TRUE.   , RealArray = U_N_2D_local)
END ASSOCIATE

DEALLOCATE(StrVarNames)
DEALLOCATE(U_N_2D_local)

GETTIME(EndT)
CALL DisplayMessageAndTime(EndT-StartT, 'DONE', DisplayDespiteLB=.TRUE., DisplayLine=.FALSE.)

END SUBROUTINE WriteElectroMagneticPICFieldToHDF5

#endif /*!(USE_FV) || (USE_HDG)*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
#endif /*defined(PARTICLES)*/
END MODULE MOD_HDF5_Output_Particles_PIC