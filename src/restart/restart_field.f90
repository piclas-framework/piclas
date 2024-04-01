!==================================================================================================================================
! Copyright (c) 2010 - 2022 Prof. Claus-Dieter Munz and Prof. Stefanos Fasoulas
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
#if USE_PETSC
#include "petsc/finclude/petsc.h"
#endif

!===================================================================================================================================
!> Module contains the routines for load balancing
!===================================================================================================================================
MODULE MOD_Restart_Field
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

!#if USE_LOADBALANCE
#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
INTERFACE FieldRestart
  MODULE PROCEDURE FieldRestart
END INTERFACE

PUBLIC :: FieldRestart
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/
!#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

CONTAINS


#if !((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))
!#if USE_LOADBALANCE
SUBROUTINE FieldRestart()
!===================================================================================================================================
! routine perfoming the field restart
!===================================================================================================================================
! USED MODULES
USE MOD_Globals
USE MOD_PreProc
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
!USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
#endif /*USE_LOADBALANCE*/
USE MOD_IO_HDF5
USE MOD_Restart_Vars       ,ONLY: N_Restart,RestartFile,RestartNullifySolution
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D,ChangeBasis2D
USE MOD_HDF5_Input         ,ONLY: OpenDataFile,CloseDataFile,ReadArray,ReadAttribute,GetDataSize
USE MOD_HDF5_Input         ,ONLY: DatasetExists
USE MOD_HDF5_Output        ,ONLY: FlushHDF5
USE MOD_Interpolation_Vars ,ONLY: N_Inter,PREF_VDM
#ifdef PP_POIS
USE MOD_Equation_Vars      ,ONLY: Phi
#elif USE_HDG
USE MOD_Interpolation_Vars ,ONLY: NMax
#else /*not PP_POIS and not USE_HDG*/
USE MOD_DG_Vars            ,ONLY: U_N
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
#endif /*PP_POIS*/
#if defined(PARTICLES)
USE MOD_Particle_Restart   ,ONLY: ParticleRestart
#endif /*defined(PARTICLES)*/
#if USE_HDG
USE MOD_Restart_Tools      ,ONLY: RecomputeLambda
USE MOD_HDG_Vars           ,ONLY: HDG_Surf_N, nGP_face
USE MOD_HDG                ,ONLY: RestartHDG
USE MOD_Mesh_Vars          ,ONLY: nSides,GlobalUniqueSideID,MortarType,SideToElem
USE MOD_StringTools        ,ONLY: set_formatting,clear_formatting
USE MOD_Mappings           ,ONLY: CGNS_SideToVol2
USE MOD_Mesh_Vars          ,ONLY: firstMortarInnerSide,lastMortarInnerSide,MortarInfo
USE MOD_Mesh_Vars          ,ONLY: lastMPISide_MINE
#if USE_MPI
USE MOD_MPI_Vars           ,ONLY: RecRequest_U,SendRequest_U
!USE MOD_MPI                ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_MPI                ,ONLY: StartReceiveMPISurfDataType,StartSendMPISurfDataType,FinishExchangeMPISurfDataType, Mask_MPIsides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_HDG                ,ONLY: SynchronizeVoltageOnEPC
USE MOD_HDG_Vars           ,ONLY: UseEPC
#if defined(PARTICLES)
USE MOD_Equation_Tools     ,ONLY: SynchronizeCPP
USE MOD_HDG                ,ONLY: SynchronizeBV
USE MOD_HDG_Vars           ,ONLY: UseBiasVoltage,UseCoupledPowerPotential
! TODO: make ElemInfo available with PARTICLES=OFF and remove this preprocessor if/else as soon as possible                                                                                                     ! TODO: make ElemInfo available with PARTICLES=OFF and remove this preprocessor if/else as soon as possible
USE MOD_Mesh_Vars          ,ONLY: SideToNonUniqueGlobalSide,nElems
USE MOD_LoadBalance_Vars   ,ONLY: MPInSideSend,MPInSideRecv,MPIoffsetSideSend,MPIoffsetSideRecv
USE MOD_Particle_Mesh_Vars ,ONLY: ElemInfo_Shared
USE MOD_HDG_Vars           ,ONLY: lambdaLB
#endif /*defined(PARTICLES)*/
#endif /*USE_LOADBALANCE*/
#if USE_PETSC
#if USE_LOADBALANCE
USE MOD_HDG                ,ONLY: SynchronizeChargeOnFPC
USE MOD_HDG_Vars           ,ONLY: UseFPC
#endif /*USE_LOADBALANCE*/
USE PETSc
USE MOD_HDG_Vars           ,ONLY: lambda_petsc,PETScGlobal,PETScLocalToSideID,nPETScUniqueSides
#endif
USE MOD_Mesh_Vars          ,ONLY: N_SurfMesh
#else /*USE_HDG*/
! Non-HDG stuff
USE MOD_Interpolation_Vars ,ONLY: Nmax
USE MOD_LoadBalance_Vars   ,ONLY: nElemsOld,offsetElemOld
USE MOD_PML_Vars           ,ONLY: DoPML,PMLToElem,nPMLElems,PMLnVar
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars   ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
#endif /*USE_LOADBALANCE*/
#endif /*USE_HDG*/
USE MOD_Mesh_Vars          ,ONLY: nElems,OffsetElem
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=IK)                   :: Nres,nVar
INTEGER(KIND=IK)                   :: OffsetElemTmp,PP_nElemsTmp
#if USE_HDG
LOGICAL                            :: DG_SolutionLambdaExists
!LOGICAL                            :: DG_SolutionExists
INTEGER                            :: SideID,iSide,MinGlobalSideID,MaxGlobalSideID,NSideMin,iVar
REAL,ALLOCATABLE                   :: ExtendedLambda(:,:,:)
INTEGER                            :: p,q,r,rr,pq(1:2)
INTEGER                            :: iLocSide,iLocSide_NB,iLocSide_master
INTEGER                            :: iMortar,MortarSideID,nMortars
REAL,ALLOCATABLE                   :: lambdaLBTmp(:,:,:)        !< lambda, ((PP_N+1)^2,nSides)
REAL,ALLOCATABLE                   :: tmp2(:,:,:),tmp3(:,:,:),lambdaloc(:,:)
#if USE_LOADBALANCE
#if defined(PARTICLES)
! TODO: make ElemInfo available with PARTICLES=OFF and remove this preprocessor if/else as soon as possible
INTEGER                            :: NonUniqueGlobalSideID
#endif /*defined(PARTICLES)*/
!INTEGER           :: checkRank
#endif /*USE_LOADBALANCE*/
#if USE_PETSC
INTEGER                            :: PETScLocalID
PetscErrorCode                     :: ierr
#endif
#else /*USE_HDG*/
INTEGER                            :: iElem
#if USE_LOADBALANCE
REAL,ALLOCATABLE                   :: UTmp(:,:,:,:,:)
#endif /*USE_LOADBALANCE*/
REAL,ALLOCATABLE                   :: U_local(:,:,:,:,:)
!REAL,ALLOCATABLE                   :: U_local2(:,:,:,:,:)
INTEGER                            :: iPML
INTEGER(KIND=IK)                   :: PMLnVarTmp
INTEGER                            :: iPMLElem
REAL,ALLOCATABLE                   :: U(:,:,:,:,:)
INTEGER                            :: i,j,k
#endif /*USE_HDG*/
#if USE_LOADBALANCE
!INTEGER,ALLOCATABLE                :: N_DG_Tmp(:)
#if defined(PARTICLES) || !(USE_HDG)
! TODO: make ElemInfo available with PARTICLES=OFF and remove this preprocessor if/else as soon as possible
! Custom data type
INTEGER                            :: MPI_LENGTH(1),MPI_TYPE(1),MPI_STRUCT
INTEGER(KIND=MPI_ADDRESS_KIND)     :: MPI_DISPLACEMENT(1)
#endif /*USE_LOADBALANCE*/
#endif /*defined(PARTICLES) || !(USE_HDG)*/
#ifdef PP_POIS
#elif USE_HDG
#else /*not PP_POIS and not USE_HDG*/
REAL,ALLOCATABLE                   :: Uloc(:,:,:,:)
INTEGER                            :: Nloc
#endif /*PP_POIS*/
!===================================================================================================================================

! ===========================================================================
! Distribute or read the field solution
! ===========================================================================

#if USE_LOADBALANCE
IF(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))THEN

    ! ! p-adaption: Re-distribute N_DG
    ! ALLOCATE(N_DG_Tmp(1:nElems))
    ! ASSOCIATE (&
    !         counts_send  => INT(MPInElemSend     ) ,&
    !         disp_send    => INT(MPIoffsetElemSend) ,&
    !         counts_recv  => INT(MPInElemRecv     ) ,&
    !         disp_recv    => INT(MPIoffsetElemRecv))
    !   MPI_LENGTH       = 1
    !   MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    !   MPI_TYPE         = MPI_INTEGER_INT_KIND
    !   CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    !   CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
    !   CALL MPI_ALLTOALLV(N_DG_Mapping(2,1+offSetElem:nElems+offSetElem),counts_send,disp_send,MPI_STRUCT,N_DG_Tmp,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
    !   CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
    ! END ASSOCIATE

#if USE_HDG
#if USE_PETSC
  ! FPC: The MPI root process distributes the information among the sub-communicator processes for each FPC
  !      (before and after load balancing, the root process is always part of each sub-communicator group)
  IF(UseFPC) CALL SynchronizeChargeOnFPC()
#endif /*USE_PETSC*/
  ! EPC: The MPI root process distributes the information among the sub-communicator processes for each EPC
  !      (before and after load balancing, the root process is always part of each sub-communicator group)
  IF(UseEPC) CALL SynchronizeVoltageOnEPC()

#if defined(PARTICLES)
  ! Coupled Bias voltage (BV): The MPI root process distributes the information among the sub-communicator processes
  !      (before and after load balancing, the root process is always part of each sub-communicator group)
  IF(UseBiasVoltage) CALL SynchronizeBV()

  ! Coupled Power Potential (CPP): The MPI root process distributes the information among the sub-communicator processes
  !      (before and after load balancing, the root process is always part of each sub-communicator group)
  IF(UseCoupledPowerPotential) CALL SynchronizeCPP()

  !checkRank=MIN(3,nProcessors-1)
  ! Store lambda solution on global non-unique array for MPI communication
  ASSOCIATE( firstSide => ElemInfo_Shared(ELEM_FIRSTSIDEIND,offsetElem+1) + 1       ,&
             lastSide  => ElemInfo_Shared(ELEM_LASTSIDEIND ,offsetElem    + nElems) )
         !IPWRITE(UNIT_StdOut,*) "firstSide:lastSide,Nbr,nSides =", firstSide,lastSide,lastSide-firstSide+1,nSides
    ALLOCATE(lambdaLBTmp(PP_nVar,nGP_face(NMax)+1,firstSide:lastSide)) ! +1 comes from the NSideMin info that is sent additionally
  END ASSOCIATE
       !CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)
       !IPWRITE(UNIT_StdOut,*) "MPInSideSend,MPIoffsetSideSend,MPInSideRecv,MPIoffsetSideRecv =", MPInSideSend,MPIoffsetSideSend,MPInSideRecv,MPIoffsetSideRecv
       !CALL MPI_BARRIER(MPI_COMM_PICLAS,iError)

  ASSOCIATE (&
          counts_send  => (INT(MPInSideSend     )) ,&
          disp_send    => (INT(MPIoffsetSideSend)) ,&
          counts_recv  => (INT(MPInSideRecv     )) ,&
          disp_recv    => (INT(MPIoffsetSideRecv)))
    ! Communicate PartInt over MPI
    MPI_LENGTH       = PP_nVar*(nGP_face(NMax)+1) ! +1 comes from the NSideMin info that is sent additionally
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
    CALL MPI_ALLTOALLV(lambdaLB,counts_send,disp_send,MPI_STRUCT,lambdaLBTmp,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  DEALLOCATE(lambdaLB)

  ! Loop over all sides and map the solution via the non-unique global side ID to the SideID (1 to nSides)
  DO iSide = 1, nSides
    IF(iSide.LE.lastMPISide_MINE)THEN
      NonUniqueGlobalSideID = SideToNonUniqueGlobalSide(1,iSide)
      iLocSide        = SideToElem(S2E_LOC_SIDE_ID    , iSide)
      iLocSide_master = SideToElem(S2E_LOC_SIDE_ID    , iSide)
      iLocSide_NB     = SideToElem(S2E_NB_LOC_SIDE_ID , iSide)

      ! Small virtual mortar master side (blue) is encountered with MortarType(1,iSide) = 0
      ! Blue (small mortar master) side writes as yellow (big mortar master) for consistency (you spin me right round baby right round)
      IF(MortarType(1,iSide).EQ.0)THEN
        ! check all my big mortar sides and find the one to which the small virtual is connected
        ! check all yellow (big mortar) sides
        Check1: DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide
          nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
          ! loop over all blue sides (small mortar master)
          DO iMortar=1,nMortars
            SideID = MortarInfo(MI_SIDEID,iMortar,MortarType(2,MortarSideID)) !small SideID
            ! check if for "21,22,23" and find the "3"
            ! iSide=21,22,23 is blue
            ! SideID=3 is yellow
            IF(iSide.EQ.SideID)THEN
              iLocSide_master = SideToElem(S2E_LOC_SIDE_ID,MortarSideID)
              IF(iLocSide_master.EQ.-1)THEN
                CALL abort(__STAMP__,'This big mortar side must be master')
              END IF !iLocSide.NE.-1
              EXIT Check1
            END IF ! iSide.EQ.SideID
          END DO !iMortar
        END DO Check1 !MortarSideID
      END IF ! MortarType(1,iSide).EQ.0

      ! Get NSideMin from last entry. Here no change basis is required as the NMax array is only filled up to the degree NSideMin
      NSideMin = INT(lambdaLBTmp(1,nGP_face(NMax)+1,NonUniqueGlobalSideID))

      ! Rotate data into correct orientation
      DO q=0,NSideMin
        DO p=0,NSideMin
          pq = CGNS_SideToVol2(NSideMin,p,q,iLocSide_master)
          r  = q    *(NSideMin+1)+p    +1
          rr = pq(2)*(NSideMin+1)+pq(1)+1
          !lambda(:,r:r,iSide) = lambdaLBTmp(:,rr:rr,NonUniqueGlobalSideID)
          HDG_Surf_N(iSide)%lambda(:,r:r) = lambdaLBTmp(:,rr:rr,NonUniqueGlobalSideID)
        END DO
      END DO !p,q
    END IF ! iSide.LE.lastMPISide_MINE
  END DO

  DEALLOCATE(lambdaLBTmp)

#if USE_MPI
  ! Exchange lambda MINE -> YOUR direction (as only the master sides have read the solution until now)
  DO iVar=1,PP_nVar
    CALL StartReceiveMPISurfDataType(RecRequest_U , 1 , 2)
    CALL StartSendMPISurfDataType(  SendRequest_U , 1 , 2, iVar) ! 2 = lambda
    CALL FinishExchangeMPISurfDataType(SendRequest_U, RecRequest_U, 1, 2, iVar) ! 2 = lambda
  END DO ! iVar=1,PP_nVar
#endif /*USE_MPI*/

#if USE_PETSC
  DO PETScLocalID=1,nPETScUniqueSides
    SideID=PETScLocalToSideID(PETScLocalID)
    PetscCallA(VecSetValuesBlocked(lambda_petsc,1,PETScGlobal(SideID),lambda(1,:,SideID),INSERT_VALUES,ierr))
  END DO
  PetscCallA(VecAssemblyBegin(lambda_petsc,ierr))
  PetscCallA(VecAssemblyEnd(lambda_petsc,ierr))
#endif

  CALL RestartHDG() ! calls PostProcessGradient for calculate the derivative, e.g., the electric field E

#else
  ! TODO: make ElemInfo available with PARTICLES=OFF and remove this preprocessor if/else as soon as possible
   CALL abort(__STAMP__,'TODO: make ElemInfo available with PARTICLES=OFF and remove this preprocessor if/else')
#endif /*defined(PARTICLES)*/

#else /*USE_HDG*/
  ! Only required for time discs where U is allocated
  ALLOCATE(U(PP_nVar,0:NMax,0:NMax,0:NMax,nElemsOld))
  DO iElem = 1, nElemsOld
    Nloc = N_DG_Mapping(2,iElem+offSetElemOld)
    IF(Nloc.EQ.Nmax)THEN
      U(:,:,:,:,iElem) = U_N(iElem)%U(:,:,:,:)
    ELSE
      U(:,:,:,:,iElem) = 0.
      DO k=0,Nloc
          DO i=0,Nloc
        DO j=0,Nloc
            U(:,i,j,k,iElem) = U_N(iElem)%U(:,i,j,k)
          END DO
        END DO
      END DO
    END IF ! Nloc.Eq.Nmax
  END DO ! iElem = 1, nElems

  ALLOCATE(UTmp(PP_nVar,0:NMax,0:NMax,0:NMax,nElems))
  ASSOCIATE (&
          counts_send  => (INT(MPInElemSend     )) ,&
          disp_send    => (INT(MPIoffsetElemSend)) ,&
          counts_recv  => (INT(MPInElemRecv     )) ,&
          disp_recv    => (INT(MPIoffsetElemRecv)))
    ! Communicate PartInt over MPI
    MPI_LENGTH       = PP_nVar*(NMax+1)*(NMax+1)*(NMax+1)
    MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
    MPI_TYPE         = MPI_DOUBLE_PRECISION
    CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
    CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
    CALL MPI_ALLTOALLV(U,counts_send,disp_send,MPI_STRUCT,UTmp,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
    CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
  END ASSOCIATE
  CALL MOVE_ALLOC(UTmp,U)

  DEALLOCATE(U_N)
  ! the local DG solution in physical and reference space
  ALLOCATE(U_N(1:nElems))
  DO iElem = 1, nElems
    Nloc = N_DG_Mapping(2,iElem+offSetElem)
    ALLOCATE(U_N(iElem)%U(PP_nVar,0:Nloc,0:Nloc,0:Nloc))
    DO k=0,Nloc
      DO j=0,Nloc
        DO i=0,Nloc
          U_N(iElem)%U(:,i,j,k) = U(:,i,j,k,iElem)
        END DO
      END DO
    END DO
    ALLOCATE(U_N(iElem)%Ut(PP_nVar,0:Nloc,0:Nloc,0:Nloc))
    U_N(iElem)%Ut = 0.
  END DO ! iElem = 1, PP_nElems

  ! Re-allocate the PML solution arrays
  IF(DoPML)THEN
    DO iPMLElem=1,nPMLElems
      iElem = PMLToElem(iPMLElem)
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      ALLOCATE(U_N(iElem)%U2(PMLnVar,0:Nloc,0:Nloc,0:Nloc))
      U_N(iElem)%U2 = 0.
      ALLOCATE(U_N(iElem)%U2t(PMLnVar,0:Nloc,0:Nloc,0:Nloc))
      U_N(iElem)%U2t = 0.
    END DO
  END IF ! DoPML
  DEALLOCATE(U)
#endif /*USE_HDG*/

  ! Update N_DG and delete N_DG_Tmp
  !CALL MOVE_ALLOC(N_DG_Tmp,N_DG)
  !DO iElem = 1, nElemsOld
  !  Nloc = N_DG_Mapping(2,iElem+offSetElem)
  !END DO
  !DEALLOCATE(N_DG_Tmp)

ELSE ! normal restart
#endif /*USE_LOADBALANCE*/

  IF(N_Restart.LT.1) CALL abort(__STAMP__,'N_Restart<1 is not allowed. Check correct initailisation of N_Restart!')

  ! Temp. vars for integer KIND=8 possibility
  Nres          = INT(N_Restart,IK)
  OffsetElemTmp = INT(OffsetElem,IK)
  nVar          = INT(PP_nVar,IK)
  PP_nElemsTmp  = INT(PP_nElems,IK)
#if !(USE_HDG)
  PMLnVarTmp    = INT(PMLnVar,IK)
#endif /*not USE_HDG*/




  ! ===========================================================================
  ! Read the field solution
  ! ===========================================================================
  IF(RestartNullifySolution)THEN ! Open the restart file and neglect the DG solution (only read particles if present)
    SWRITE(UNIT_stdOut,*)'Restarting from File: ',TRIM(RestartFile),' (but without reading the DG solution)'
    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
  ELSE ! Use the solution in the restart file
    SWRITE(UNIT_stdOut,*)'Restarting from File: ',TRIM(RestartFile)
    CALL OpenDataFile(RestartFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_PICLAS)
    ! Read in time from restart file
    !CALL ReadAttribute(File_ID,'Time',1,RealScalar=RestartTime)
    ! Read in state


    SWRITE(UNIT_stdOut,*)'Interpolating solution from restart grid with N=',N_restart,' to computational grid with N=',PP_N



#ifdef PP_POIS
#if (PP_nVar==8)
    CALL ReadArray('DG_SolutionE',5,(/PP_nVar,Nres+1_IK,Nres+1_IK,Nres+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
    CALL ReadArray('DG_SolutionPhi',5,(/4_IK,Nres+1_IK,Nres+1_IK,Nres+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=Phi)
#else
    CALL ReadArray('DG_SolutionE',5,(/PP_nVar,Nres+1_IK,Nres+1_IK,Nres+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
    CALL ReadArray('DG_SolutionPhi',5,(/PP_nVar,Nres+1_IK,Nres+1_IK,Nres+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=Phi)
#endif /*PP_nVar==8*/
#elif USE_HDG
    ! TODO: Do we need this for the HDG solver?
    !CALL DatasetExists(File_ID,'DG_Solution',DG_SolutionExists)
    !IF(DG_SolutionExists)THEN
    !  ALLOCATE(U(PP_nVar,0:Nres,0:Nres,0:Nres,PP_nElemsTmp))
    !  CALL ReadArray('DG_Solution' ,5,(/PP_nVar,Nres+1_IK,Nres+1_IK,Nres+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
    !ELSE
    !  CALL abort(__STAMP__,'DG_Solution does not exist')
    !END IF

    ! Read HDG lambda solution (sorted in ascending global unique side ID ordering)
    CALL DatasetExists(File_ID,'DG_SolutionLambda',DG_SolutionLambdaExists)

    IF(DG_SolutionLambdaExists)THEN
      MinGlobalSideID = HUGE(1)
      MaxGlobalSideID = -1
      DO iSide = 1, nSides
        MaxGlobalSideID = MERGE(ABS(GlobalUniqueSideID(iSide)) , MaxGlobalSideID , ABS(GlobalUniqueSideID(iSide)).GT.MaxGlobalSideID)
        MinGlobalSideID = MERGE(ABS(GlobalUniqueSideID(iSide)) , MinGlobalSideID , ABS(GlobalUniqueSideID(iSide)).LT.MinGlobalSideID)
      END DO

        ASSOCIATE( &
              ExtendedOffsetSide => INT(MinGlobalSideID-1,IK)                 ,&
              ExtendednSides     => INT(MaxGlobalSideID-MinGlobalSideID+1,IK) ,&
              nGP_faceNMax       => INT(nGP_face(NMax),IK)                    ,&
              PP_nVarTmp         => INT(PP_nVar,IK))
          !ALLOCATE(ExtendedLambda(PP_nVar,nGP_face,MinGlobalSideID:MaxGlobalSideID))
          ALLOCATE(ExtendedLambda(PP_nVar,nGP_faceNMax,1:ExtendednSides))
          ExtendedLambda = HUGE(1.)
          !lambda = HUGE(1.)
          ! WARGNING: Data read in overlapping mode, therefore ReadNonOverlap_opt=F
          CALL ReadArray('DG_SolutionLambda',3,(/PP_nVarTmp,nGP_faceNMax,ExtendednSides/),ExtendedOffsetSide,3,RealArray=ExtendedLambda&
#if USE_MPI
              ,ReadNonOverlap_opt=.FALSE.&
#endif /*USE_MPI*/
              )

          ! Allocate temporary arrays
          ALLOCATE(tmp2(1:PP_nVar,0:Nres,0:Nres))
          ALLOCATE(tmp3(1:PP_nVar,0:Nres,0:Nres))
          ALLOCATE(lambdaloc(1:PP_nVar,1:nGP_face(Nres)))

          DO iSide = 1, nSides
            IF(iSide.LE.lastMPISide_MINE)THEN
              iLocSide        = SideToElem(S2E_LOC_SIDE_ID    , iSide)
              iLocSide_master = SideToElem(S2E_LOC_SIDE_ID    , iSide)
              iLocSide_NB     = SideToElem(S2E_NB_LOC_SIDE_ID , iSide)

              ! Small virtual mortar master side (blue) is encountered with MortarType(1,iSide) = 0
              ! Blue (small mortar master) side writes as yellow (big mortar master) for consistency (you spin me right round baby right round)
              IF(MortarType(1,iSide).EQ.0)THEN
                ! check all my big mortar sides and find the one to which the small virtual is connected
                ! check all yellow (big mortar) sides
                Check1: DO MortarSideID=firstMortarInnerSide,lastMortarInnerSide
                  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
                  ! loop over all blue sides (small mortar master)
                  DO iMortar=1,nMortars
                    SideID = MortarInfo(MI_SIDEID,iMortar,MortarType(2,MortarSideID)) !small SideID
                    ! check if for "21,22,23" and find the "3"
                    ! iSide=21,22,23 is blue
                    ! SideID=3 is yellow
                    IF(iSide.EQ.SideID)THEN
                      iLocSide_master = SideToElem(S2E_LOC_SIDE_ID,MortarSideID)
                      IF(iLocSide_master.EQ.-1)THEN
                        CALL abort(__STAMP__,'This big mortar side must be master')
                      END IF !iLocSide.NE.-1
                      EXIT Check1
                    END IF ! iSide.EQ.SideID
                  END DO !iMortar
                END DO Check1 !MortarSideID
              END IF ! MortarType(1,iSide).EQ.0

              ! Read lambda from h5 on Nres
              DO q=0,Nres
                DO p=0,Nres
                  pq = CGNS_SideToVol2(Nres,p,q,iLocSide_master)
                  r  = q    *(Nres+1)+p    +1
                  rr = pq(2)*(Nres+1)+pq(1)+1
                  lambdaloc(:,r:r) = ExtendedLambda(:,rr:rr,GlobalUniqueSideID(iSide)-ExtendedOffsetSide)
                END DO
              END DO !p,q

              ! Map from 1D to 2D array
              DO iVar=1,PP_nVar
                tmp2(iVar:iVar,0:NSideMin,0:NSideMin) = RESHAPE(lambdaloc(iVar,:),(/1,NSideMin+1,NSideMin+1/))
              END DO ! iVar=1,PP_nVar

              ! Get NSideMin
              NSideMin = N_SurfMesh(iSide)%NSideMin

              ! Map lambda from Nres to NSideMin
              IF(NSideMin.EQ.N_Restart)THEN ! N is equal
                tmp3(1:PP_nVar,0:NSideMin,0:NSideMin) = tmp2(1:PP_nVar,0:Nres,0:Nres)
              ELSEIF(NSideMin.GT.N_Restart)THEN ! N increases: Simply interpolate the lower polynomial degree solution
                CALL ChangeBasis2D(PP_nVar, N_Restart, NSideMin, PREF_VDM(N_Restart, NSideMin)%Vdm, tmp2(1:PP_nVar,0:Nres,0:Nres), tmp3(1:PP_nVar,0:NSideMin,0:NSideMin))
              ELSE ! N reduces: This requires an intermediate modal basis
                ! Switch to Legendre basis
                CALL ChangeBasis2D(PP_nVar, N_Restart, N_Restart, N_Inter(N_Restart)%sVdm_Leg, tmp2(1:PP_nVar,0:Nres,0:Nres), tmp2(1:PP_nVar,0:Nres,0:Nres))
                ! Switch back to nodal basis but cut-off the higher-order DOFs
                CALL ChangeBasis2D(PP_nVar, NSideMin, NSideMin, N_Inter(NSideMin)%Vdm_Leg, tmp2(1:PP_nVar,0:NSideMin,0:NSideMin), tmp3(1:PP_nVar,0:NSideMin,0:NSideMin))
              END IF ! NSideMin.EQ.N_Restart

              ! Map back from 2D to 1D array
              DO iVar=1,PP_nVar
                HDG_Surf_N(iSide)%lambda(iVar,:) = RESHAPE(tmp3(1,0:NSideMin,0:NSideMin),(/nGP_face(NSideMin)/))
              END DO ! iVar=1,PP_nVar

            END IF ! iSide.LE.lastMPISide_MINE
          END DO
          DEALLOCATE(ExtendedLambda)
          DEALLOCATE(tmp2)
          DEALLOCATE(tmp3)
          DEALLOCATE(lambdaloc)
        END ASSOCIATE

#if USE_MPI
      ! Exchange lambda MINE -> YOUR direction (as only the master sides have read the solution until now)
      DO iVar=1,PP_nVar
        CALL StartReceiveMPISurfDataType(RecRequest_U , 1 , 2)
        CALL StartSendMPISurfDataType(  SendRequest_U , 1 , 2, iVar) ! 2 = lambda
        CALL FinishExchangeMPISurfDataType(SendRequest_U, RecRequest_U, 1, 2, iVar) ! 2 = lambda
      END DO ! iVar=1,PP_nVar
#endif /*USE_MPI*/

#if USE_PETSC
        DO PETScLocalID=1,nPETScUniqueSides
          SideID=PETScLocalToSideID(PETScLocalID)
          PetscCallA(VecSetValuesBlocked(lambda_petsc,1,PETScGlobal(SideID),lambda(1,:,SideID),INSERT_VALUES,ierr))
        END DO
        PetscCallA(VecAssemblyBegin(lambda_petsc,ierr))
        PetscCallA(VecAssemblyEnd(lambda_petsc,ierr))
#endif

        CALL RestartHDG() ! calls PostProcessGradient for calculate the derivative, e.g., the electric field E

    ELSE
      DO SideID = 1, nSides
        HDG_Surf_N(SideID)%lambda=0.
      END DO ! SideID = 1, nSides
    END IF

#else /*not PP_POIS and not USE_HDG*/
    ALLOCATE(U(1:nVar,0:Nres,0:Nres,0:Nres,PP_nElemsTmp))
    ALLOCATE(Uloc(1:nVar,0:Nres,0:Nres,0:Nres))
    CALL ReadArray('DG_Solution',5,(/nVar,Nres+1_IK,Nres+1_IK,Nres+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
    DO iElem = 1, nElems
      Nloc = N_DG_Mapping(2,iElem+offSetElem)
      IF(Nloc.EQ.N_Restart)THEN ! N is equal
        U_N(iElem)%U(1:nVar,0:Nres,0:Nres,0:Nres) = U(1:nVar,0:Nres,0:Nres,0:Nres,iElem)
      ELSEIF(Nloc.GT.N_Restart)THEN ! N increases
        CALL ChangeBasis3D(PP_nVar, N_Restart, Nloc, PREF_VDM(N_Restart, Nloc)%Vdm, U(1:nVar,0:Nres,0:Nres,0:Nres,iElem), U_N(iElem)%U(1:nVar,0:Nloc,0:Nloc,0:Nloc))
      ELSE ! N reduces
        !transform the slave side to the same degree as the master: switch to Legendre basis
        CALL ChangeBasis3D(PP_nVar, N_Restart, N_Restart, N_Inter(N_Restart)%sVdm_Leg, U(1:nVar,0:Nres,0:Nres,0:Nres,iElem), Uloc)
        ! switch back to nodal basis
        CALL ChangeBasis3D(PP_nVar, Nloc, Nloc, N_Inter(Nloc)%Vdm_Leg, Uloc(1:nVar,0:Nloc,0:Nloc,0:Nloc), U_N(iElem)%U(1:nVar,0:Nloc,0:Nloc,0:Nloc))
      END IF ! Nloc.EQ.N_Restart
    END DO ! iElem = 1, nElems
    DEALLOCATE(Uloc)
    IF(DoPML)THEN
      ALLOCATE(U_local(PMLnVar,0:Nres,0:Nres,0:Nres,nElems))
      ALLOCATE(Uloc(1:PMLnVar,0:Nres,0:Nres,0:Nres))
      CALL ReadArray('PML_Solution',5,(/INT(PMLnVar,IK),Nres+1_IK,Nres+1_IK,Nres+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iPML=1,nPMLElems
        iElem = PMLToElem(iPML)
        Nloc  = N_DG_Mapping(2,iElem+offSetElem)
        IF(Nloc.EQ.N_Restart)THEN ! N is equal
          U_N(iElem)%U2(1:PMLnVar,0:Nres,0:Nres,0:Nres) = U_local(1:PMLnVar,0:Nres,0:Nres,0:Nres,iElem)
        ELSEIF(Nloc.GT.N_Restart)THEN ! N increases
          CALL ChangeBasis3D(PP_nVar, N_Restart, Nloc, PREF_VDM(N_Restart, Nloc)%Vdm, U_local(1:PMLnVar,0:Nres,0:Nres,0:Nres,iElem), U_N(iElem)%U2(1:PMLnVar,0:Nloc,0:Nloc,0:Nloc))
        ELSE ! N reduces
          !transform the slave side to the same degree as the master: switch to Legendre basis
          CALL ChangeBasis3D(PP_nVar, N_Restart, N_Restart, N_Inter(N_Restart)%sVdm_Leg, U_local(1:PMLnVar,0:Nres,0:Nres,0:Nres,iElem), Uloc)
          ! switch back to nodal basis
          CALL ChangeBasis3D(PP_nVar, Nloc, Nloc, N_Inter(Nloc)%Vdm_Leg, Uloc(1:PMLnVar,0:Nloc,0:Nloc,0:Nloc), U_N(iElem)%U2(1:PMLnVar,0:Nloc,0:Nloc,0:Nloc))
        END IF ! Nloc.EQ.N_Restart
      END DO ! iPML
      DEALLOCATE(U_local)
      DEALLOCATE(Uloc)
    END IF ! DoPML
#endif /*PP_POIS*/
    !CALL ReadState(RestartFile,nVar,PP_N,PP_nElems,U)









  END IF ! IF(.NOT. RestartNullifySolution)
  CALL CloseDataFile()

#if USE_LOADBALANCE
END IF ! PerformLoadBalance
#endif /*USE_LOADBALANCE*/

END SUBROUTINE FieldRestart
!#endif /*USE_LOADBALANCE*/
#endif /*!((PP_TimeDiscMethod==4) || (PP_TimeDiscMethod==300) || (PP_TimeDiscMethod==400))*/

END MODULE MOD_Restart_Field
