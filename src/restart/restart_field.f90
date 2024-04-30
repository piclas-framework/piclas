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
INTERFACE FieldRestart
  MODULE PROCEDURE FieldRestart
END INTERFACE

PUBLIC :: FieldRestart
!#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

CONTAINS

!#if USE_LOADBALANCE
SUBROUTINE FieldRestart()
!===================================================================================================================================
! routine perfoming the field restart
!===================================================================================================================================
! USED MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars                ,ONLY: U
#if USE_LOADBALANCE
USE MOD_LoadBalance_Vars       ,ONLY: PerformLoadBalance,UseH5IOLoadBalance
#endif /*USE_LOADBALANCE*/
USE MOD_IO_HDF5
USE MOD_Restart_Vars           ,ONLY: N_Restart,RestartFile,InterpolateSolution,RestartNullifySolution
USE MOD_ChangeBasis            ,ONLY: ChangeBasis3D
USE MOD_HDF5_Input             ,ONLY: OpenDataFile,CloseDataFile,ReadArray,ReadAttribute,GetDataSize
USE MOD_HDF5_Input             ,ONLY: DatasetExists
USE MOD_HDF5_Output            ,ONLY: FlushHDF5
#ifdef PP_POIS
USE MOD_Equation_Vars          ,ONLY: Phi
#endif /*PP_POIS*/
#if defined(PARTICLES)
USE MOD_Particle_Restart       ,ONLY: ParticleRestart
#endif /*defined(PARTICLES)*/
#if USE_HDG
USE MOD_Restart_Tools          ,ONLY: RecomputeLambda
USE MOD_HDG_Vars               ,ONLY: lambda, nGP_face
USE MOD_HDG                    ,ONLY: RestartHDG
USE MOD_Mesh_Vars              ,ONLY: nSides,GlobalUniqueSideID,MortarType,SideToElem
USE MOD_StringTools            ,ONLY: set_formatting,clear_formatting
USE MOD_Mappings               ,ONLY: CGNS_SideToVol2
USE MOD_Mesh_Vars              ,ONLY: firstMortarInnerSide,lastMortarInnerSide,MortarInfo
USE MOD_Mesh_Vars              ,ONLY: lastMPISide_MINE
#if USE_MPI
USE MOD_MPI_Vars               ,ONLY: RecRequest_U,SendRequest_U
USE MOD_MPI                    ,ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_HDG                    ,ONLY: SynchronizeVoltageOnEPC
USE MOD_HDG_Vars               ,ONLY: UseEPC
#if defined(PARTICLES)
USE MOD_Equation_Tools         ,ONLY: SynchronizeCPP
USE MOD_HDG                    ,ONLY: SynchronizeBV
USE MOD_HDG_Vars               ,ONLY: UseBiasVoltage,UseCoupledPowerPotential
! TODO: make ElemInfo available with PARTICLES=OFF and remove this preprocessor if/else as soon as possible
USE MOD_Mesh_Vars              ,ONLY: SideToNonUniqueGlobalSide,nElems
USE MOD_LoadBalance_Vars       ,ONLY: MPInSideSend,MPInSideRecv,MPIoffsetSideSend,MPIoffsetSideRecv
USE MOD_Particle_Mesh_Vars     ,ONLY: ElemInfo_Shared
USE MOD_HDG_Vars               ,ONLY: lambdaLB
#endif /*defined(PARTICLES)*/
#endif /*USE_LOADBALANCE*/
#if USE_PETSC
#if USE_LOADBALANCE
USE MOD_HDG                    ,ONLY: SynchronizeChargeOnFPC
USE MOD_HDG_Vars               ,ONLY: UseFPC
#endif /*USE_LOADBALANCE*/
USE PETSc
USE MOD_HDG_Vars               ,ONLY: lambda_petsc,PETScGlobal,PETScLocalToSideID,nPETScUniqueSides
#endif
#else /*USE_HDG*/
! Non-HDG stuff
USE MOD_PML_Vars               ,ONLY: DoPML,PMLToElem,U2,nPMLElems,PMLnVar
USE MOD_Restart_Vars           ,ONLY: Vdm_GaussNRestart_GaussN
#if USE_LOADBALANCE
USE MOD_Mesh_Vars              ,ONLY: nElems
USE MOD_LoadBalance_Vars       ,ONLY: MPInElemSend,MPInElemRecv,MPIoffsetElemSend,MPIoffsetElemRecv
#endif /*USE_LOADBALANCE*/
#endif /*USE_HDG*/
USE MOD_Mesh_Vars              ,ONLY: OffsetElem
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=IK)                   :: PP_NTmp,OffsetElemTmp,PP_nVarTmp,PP_nElemsTmp,N_RestartTmp
#if USE_HDG
LOGICAL                            :: DG_SolutionLambdaExists
LOGICAL                            :: DG_SolutionUExists
INTEGER                            :: SideID,iSide,MinGlobalSideID,MaxGlobalSideID
REAL,ALLOCATABLE                   :: ExtendedLambda(:,:,:)
INTEGER                            :: p,q,r,rr,pq(1:2)
INTEGER                            :: iLocSide,iLocSide_NB,iLocSide_master
INTEGER                            :: iMortar,MortarSideID,nMortars
#if USE_LOADBALANCE
#if defined(PARTICLES)
! TODO: make ElemInfo available with PARTICLES=OFF and remove this preprocessor if/else as soon as possible
REAL,ALLOCATABLE                   :: lambdaLBTmp(:,:,:)        !< lambda, ((PP_N+1)^2,nSides)
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
REAL,ALLOCATABLE                   :: U_local2(:,:,:,:,:)
INTEGER                            :: iPML
INTEGER(KIND=IK)                   :: PMLnVarTmp
#endif /*USE_HDG*/
#if USE_LOADBALANCE
#if defined(PARTICLES) || !(USE_HDG)
! TODO: make ElemInfo available with PARTICLES=OFF and remove this preprocessor if/else as soon as possible
! Custom data type
INTEGER                            :: MPI_LENGTH(1),MPI_TYPE(1),MPI_STRUCT
INTEGER(KIND=MPI_ADDRESS_KIND)     :: MPI_DISPLACEMENT(1)
#endif /*defined(PARTICLES) || !(USE_HDG)*/
#endif /*USE_LOADBALANCE*/
!===================================================================================================================================

! ===========================================================================
! Distribute or read the field solution
! ===========================================================================

#if USE_LOADBALANCE
IF(PerformLoadBalance.AND.(.NOT.UseH5IOLoadBalance))THEN
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
    ALLOCATE(lambdaLBTmp(PP_nVar,nGP_face,firstSide:lastSide))
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
    MPI_LENGTH       = PP_nVar*nGP_face
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

      ! Rotate data into correct orientation
      DO q=0,PP_N
        DO p=0,PP_N
          pq = CGNS_SideToVol2(PP_N,p,q,iLocSide_master)
          r  = q    *(PP_N+1)+p    +1
          rr = pq(2)*(PP_N+1)+pq(1)+1
          lambda(:,r:r,iSide) = lambdaLBTmp(:,rr:rr,NonUniqueGlobalSideID)
        END DO
      END DO !p,q
    END IF ! iSide.LE.lastMPISide_MINE
  END DO

  DEALLOCATE(lambdaLBTmp)

#if USE_MPI
  ! Exchange lambda MINE -> YOUR direction (as only the master sides have read the solution until now)
  CALL StartReceiveMPIData(1,lambda,1,nSides, RecRequest_U,SendID=1) ! Receive YOUR
  CALL StartSendMPIData(   1,lambda,1,nSides,SendRequest_U,SendID=1) ! Send MINE
  CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=1)
#endif /*USE_MPI*/

#if USE_PETSC
  DO PETScLocalID=1,nPETScUniqueSides
    SideID=PETScLocalToSideID(PETScLocalID)
    PetscCallA(VecSetValuesBlocked(lambda_petsc,1,PETScGlobal(SideID),lambda(1,:,SideID),INSERT_VALUES,ierr))
  END DO
  PetscCallA(VecAssemblyBegin(lambda_petsc,ierr))
  PetscCallA(VecAssemblyEnd(lambda_petsc,ierr))
#endif

  CALL RestartHDG(U) ! calls PostProcessGradient for calculate the derivative, e.g., the electric field E

#else
  ! TODO: make ElemInfo available with PARTICLES=OFF and remove this preprocessor if/else as soon as possible
  lambda=0.
#endif /*defined(PARTICLES)*/

#else /*USE_HDG*/
  ! Only required for time discs where U is allocated
  IF(ALLOCATED(U))THEN
    ALLOCATE(UTmp(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems))
    ASSOCIATE (&
            counts_send  => (INT(MPInElemSend     )) ,&
            disp_send    => (INT(MPIoffsetElemSend)) ,&
            counts_recv  => (INT(MPInElemRecv     )) ,&
            disp_recv    => (INT(MPIoffsetElemRecv)))
      ! Communicate PartInt over MPI
      MPI_LENGTH       = PP_nVar*(PP_N+1)*(PP_N+1)*(PP_N+1)
      MPI_DISPLACEMENT = 0  ! 0*SIZEOF(MPI_SIZE)
      MPI_TYPE         = MPI_DOUBLE_PRECISION
      CALL MPI_TYPE_CREATE_STRUCT(1,MPI_LENGTH,MPI_DISPLACEMENT,MPI_TYPE,MPI_STRUCT,iError)
      CALL MPI_TYPE_COMMIT(MPI_STRUCT,iError)
      CALL MPI_ALLTOALLV(U,counts_send,disp_send,MPI_STRUCT,UTmp,counts_recv,disp_recv,MPI_STRUCT,MPI_COMM_PICLAS,iError)
      CALL MPI_TYPE_FREE(MPI_STRUCT,iError)
    END ASSOCIATE
    CALL MOVE_ALLOC(UTmp,U)
  END IF ! ALLOCATED(U)
#endif /*USE_HDG*/

ELSE ! normal restart
#endif /*USE_LOADBALANCE*/

  ! Temp. vars for integer KIND=8 possibility
  PP_NTmp       = INT(PP_N,IK)
  OffsetElemTmp = INT(OffsetElem,IK)
  PP_nVarTmp    = INT(PP_nVar,IK)
  PP_nElemsTmp  = INT(PP_nElems,IK)
  N_RestartTmp  = INT(N_Restart,IK)
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
    IF(.NOT. InterpolateSolution)THEN! No interpolation needed, read solution directly from file
#ifdef PP_POIS
#if (PP_nVar==8)
      CALL ReadArray('DG_SolutionE',5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
      CALL ReadArray('DG_SolutionPhi',5,(/4_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=Phi)
#else
      CALL ReadArray('DG_SolutionE',5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
      CALL ReadArray('DG_SolutionPhi',5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=Phi)
#endif
#elif USE_HDG
      CALL DatasetExists(File_ID,'DG_SolutionU',DG_SolutionUExists)
      IF(DG_SolutionUExists)THEN
        CALL ReadArray('DG_SolutionU',5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
      ELSE
        CALL ReadArray('DG_Solution' ,5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
      END IF

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
              nGP_face           => INT(nGP_face,IK)                           )
          !ALLOCATE(ExtendedLambda(PP_nVar,nGP_face,MinGlobalSideID:MaxGlobalSideID))
          ALLOCATE(ExtendedLambda(PP_nVar,nGP_face,1:ExtendednSides))
          ExtendedLambda = HUGE(1.)
          lambda = HUGE(1.)
          ! WARGNING: Data read in overlapping mode, therefore ReadNonOverlap_opt=F
          CALL ReadArray('DG_SolutionLambda',3,(/PP_nVarTmp,nGP_face,ExtendednSides/),ExtendedOffsetSide,3,RealArray=ExtendedLambda&
#if USE_MPI
              ,ReadNonOverlap_opt=.FALSE.&
#endif /*USE_MPI*/
              )

          ! Loop over all sides and map the solution via the non-unique global side ID to the SideID (1 to nSides)
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

              DO q=0,PP_N
                DO p=0,PP_N
                  pq = CGNS_SideToVol2(PP_N,p,q,iLocSide_master)
                  r  = q    *(PP_N+1)+p    +1
                  rr = pq(2)*(PP_N+1)+pq(1)+1
                  lambda(:,r:r,iSide) = ExtendedLambda(:,rr:rr,GlobalUniqueSideID(iSide)-ExtendedOffsetSide)
                END DO
              END DO !p,q
              !IPWRITE(UNIT_StdOut,*) "iSide,SUM(lambda(:,r:r,iSide)) =", iSide,SUM(lambda(:,r:r,iSide))
            END IF ! iSide.LE.lastMPISide_MINE
          END DO
          DEALLOCATE(ExtendedLambda)
        END ASSOCIATE

#if USE_MPI
        ! Exchange lambda MINE -> YOUR direction (as only the master sides have read the solution until now)
        CALL StartReceiveMPIData(1,lambda,1,nSides, RecRequest_U,SendID=1) ! Receive YOUR
        CALL StartSendMPIData(   1,lambda,1,nSides,SendRequest_U,SendID=1) ! Send MINE
        CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=1)
#endif /*USE_MPI*/

#if USE_PETSC
        DO PETScLocalID=1,nPETScUniqueSides
          SideID=PETScLocalToSideID(PETScLocalID)
          PetscCallA(VecSetValuesBlocked(lambda_petsc,1,PETScGlobal(SideID),lambda(1,:,SideID),INSERT_VALUES,ierr))
        END DO
        PetscCallA(VecAssemblyBegin(lambda_petsc,ierr))
        PetscCallA(VecAssemblyEnd(lambda_petsc,ierr))
#endif

        CALL RestartHDG(U) ! calls PostProcessGradient for calculate the derivative, e.g., the electric field E

      ELSE
        lambda=0.
      END IF

#else
      CALL ReadArray('DG_Solution',5,(/PP_nVarTmp,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),OffsetElemTmp,5,RealArray=U)
      IF(DoPML)THEN
        ALLOCATE(U_local(PMLnVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
        CALL ReadArray('PML_Solution',5,(/INT(PMLnVar,IK),PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),&
            OffsetElemTmp,5,RealArray=U_local)
        DO iPML=1,nPMLElems
          U2(:,:,:,:,iPML) = U_local(:,:,:,:,PMLToElem(iPML))
        END DO ! iPML
        DEALLOCATE(U_local)
      END IF ! DoPML
#endif
      !CALL ReadState(RestartFile,PP_nVar,PP_N,PP_nElems,U)
    ELSE! We need to interpolate the solution to the new computational grid
      SWRITE(UNIT_stdOut,*)'Interpolating solution from restart grid with N=',N_restart,' to computational grid with N=',PP_N
#ifdef PP_POIS
#if (PP_nVar==8)
      ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      CALL ReadArray('DG_SolutionE',5,(/PP_nVarTmp,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
      END DO
      DEALLOCATE(U_local)

      ALLOCATE(U_local(4,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      CALL ReadArray('DG_SolutionPhi',5,(/4_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(4,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),Phi(:,:,:,:,iElem))
      END DO
      DEALLOCATE(U_local)
#else
      ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      CALL ReadArray('DG_SolutionE',5,(/PP_nVarTmp,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
      END DO
      CALL ReadArray('DG_SolutionPhi',5,(/PP_nVarTmp,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),Phi(:,:,:,:,iElem))
      END DO
      DEALLOCATE(U_local)
#endif
#elif USE_HDG
      lambda=0.
      !CALL abort(&
          !__STAMP__&
          !,' Restart with changed polynomial degree not implemented for HDG!')
      !    ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      !    CALL ReadArray('DG_SolutionLambda',5,(/PP_nVar,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),OffsetElem,5,RealArray=U_local)
      !    DO iElem=1,PP_nElems
      !      CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
      !    END DO
      !    DEALLOCATE(U_local)
      !CALL RestartHDG(U)
#else
      ALLOCATE(U_local(PP_nVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
      CALL ReadArray('DG_Solution',5,(/PP_nVarTmp,N_RestartTmp+1_IK,N_RestartTmp+1_IK,N_RestartTmp+1_IK,PP_nElemsTmp/),&
          OffsetElemTmp,5,RealArray=U_local)
      DO iElem=1,PP_nElems
        CALL ChangeBasis3D(PP_nVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U(:,:,:,:,iElem))
      END DO
      DEALLOCATE(U_local)
      IF(DoPML)THEN
        ALLOCATE(U_local(PMLnVar,0:N_Restart,0:N_Restart,0:N_Restart,PP_nElems))
        ALLOCATE(U_local2(PMLnVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
        CALL ReadArray('PML_Solution',5,(/INT(PMLnVar,IK),PP_NTmp+1_IK,PP_NTmp+1_IK,PP_NTmp+1_IK,PP_nElemsTmp/),&
            OffsetElemTmp,5,RealArray=U_local)
        DO iElem=1,PP_nElems
          CALL ChangeBasis3D(PMLnVar,N_Restart,PP_N,Vdm_GaussNRestart_GaussN,U_local(:,:,:,:,iElem),U_local2(:,:,:,:,iElem))
        END DO
        DO iPML=1,nPMLElems
          U2(:,:,:,:,iPML) = U_local2(:,:,:,:,PMLToElem(iPML))
        END DO ! iPML
        DEALLOCATE(U_local,U_local2)
      END IF ! DoPML
#endif
      SWRITE(UNIT_stdOut,*)' DONE!'
    END IF ! IF(.NOT. InterpolateSolution)
  END IF ! IF(.NOT. RestartNullifySolution)
  CALL CloseDataFile()

#if USE_LOADBALANCE
END IF ! PerformLoadBalance
#endif /*USE_LOADBALANCE*/

END SUBROUTINE FieldRestart
!#endif /*USE_LOADBALANCE*/

END MODULE MOD_Restart_Field
