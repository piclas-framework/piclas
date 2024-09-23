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
#if USE_PETSC
#include "petsc/finclude/petsc.h"
#endif

!===================================================================================================================================
!> Module for the HDG method
!===================================================================================================================================
MODULE MOD_HDG_Tools
! MODULES
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
#if USE_HDG
PUBLIC :: CG_solver
PUBLIC :: DisplayConvergence
PUBLIC :: MatVec
PUBLIC :: evalresidual
PUBLIC :: VectorDotProductRR
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS

#if USE_HDG
!===================================================================================================================================
!> Conjugate Gradient solver
!===================================================================================================================================
SUBROUTINE CG_solver(iVar)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_HDG_Vars           ,ONLY: HDGDisplayConvergence,iteration
USE MOD_HDG_Vars           ,ONLY: EpsCG,MaxIterCG,PrecondType,useRelativeAbortCrit,OutIterCG,HDG_Surf_N
USE MOD_TimeDisc_Vars      ,ONLY: iter,IterDisplayStep
USE MOD_Mesh_Vars          ,ONLY: nSides
#if USE_MPI
!USE MOD_Mesh_Vars          ,ONLY: nMPISides_YOUR
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars           ,ONLY: MPIW8TimeField,MPIW8CountField
#endif /*defined(MEASURE_MPI_WAIT)*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(IN),OPTIONAL::iVar
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL,DIMENSION(nGP_face*nSides) :: V,Z,R
INTEGER                         :: SideID
REAL                            :: AbortCrit2
REAL                            :: omega,rr,vz,rz1,rz2,Norm_r2
REAL                            :: timestartCG,timeEndCG
!INTEGER                         :: VecSize
LOGICAL                         :: converged
#if USE_LOADBALANCE
REAL                            :: tLBStart
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)   :: CounterStart,CounterEnd
REAL(KIND=8)      :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
IF(HDGDisplayConvergence.AND.(MOD(iter,IterDisplayStep).EQ.0)) THEN
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(*,*)'CG solver start'
END IF
TimeStartCG=PICLASTIME()
#if USE_MPI
! not use MPI_YOUR sides for vector_dot_product!!!
!CALL abort(__STAMP__,'not implemented')
!VecSize=(nSides-nMPIsides_YOUR)*nGP_face(PP_N)
#else
!VecSize=nSides*nGP_face
#endif /*USE_MPI*/
IF(PRESENT(iVar)) THEN
  CALL EvalResidual(iVar)
ELSE
  CALL EvalResidual(1)
END IF

CALL VectorDotProductRR(Norm_R2) !Z=V

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

IF(useRelativeAbortCrit)THEN
#if USE_MPI
  IF(MPIroot) converged=(Norm_R2.LT.1e-16)
  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_PICLAS,iError)
#else
  converged=(Norm_R2.LT.1e-16)
#endif /*USE_MPI*/
ELSE
#if USE_MPI
  IF(MPIroot) converged=(Norm_R2.LT.EpsCG**2)
  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_PICLAS,iError)
#else
  converged=(Norm_R2.LT.EpsCG**2)
#endif /*USE_MPI*/
END IF

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(3)  = MPIW8TimeField(3) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(3) = MPIW8CountField(3) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

IF(converged) THEN !converged
!  SWRITE(*,*)'CG not needed, residual already =0'
!  SWRITE(UNIT_StdOut,'(132("-"))')
    TimeEndCG=PICLASTIME()
    iteration = 0
    IF(MPIroot) CALL DisplayConvergence(TimeEndCG-TimeStartCG, iteration, SQRT(Norm_R2))
  RETURN
END IF !converged
AbortCrit2=EpsCG**2
IF(useRelativeAbortCrit) AbortCrit2=Norm_R2*EpsCG**2

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
IF(PrecondType.NE.0) THEN
  CALL ApplyPrecond(iVar,.TRUE.) ! ApplyPrecond(R,V)
ELSE
  DO SideID = 1, nSides
    HDG_Surf_N(SideID)%V(iVar,:) = HDG_Surf_N(SideID)%R(iVar,:)
  END DO ! SideID = 1, nSides
END IF
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/
CALL VectorDotProductRV(rz1) !Z=V

! Conjugate Gradient
!IF(MPIroot) print*, '!!!!!!!!!!!!!!!!!!!!!!'
!IF(MPIroot) print*, iVar
DO iteration=1,MaxIterCG
  ! matrix vector
  IF(PRESENT(iVar)) THEN
    CALL MatVec(iVar,.TRUE.) ! CALL MatVec(V,Z, iVar)
  ELSE
    CALL MatVec(1,.TRUE.) ! CALL MatVec(V,Z)
  END IF

  CALL VectorDotProductVZ(vz)

  IF(ABS(vz).LE.0.0) CALL abort(__STAMP__,'vz is zero!')
  omega=rz1/vz

  ! Use iVar here
  DO SideID = 1, nSides
    HDG_Surf_N(SideID)%lambda(iVar,:) = HDG_Surf_N(SideID)%lambda(iVar,:) + omega * HDG_Surf_N(SideID)%V(iVar,:)
    HDG_Surf_N(SideID)%R(iVar,:)      = HDG_Surf_N(SideID)%R(iVar,:)      - omega * HDG_Surf_N(SideID)%Z(iVar,:)
  END DO ! SideID = 1, nSides
  CALL VectorDotProductRR(rr)
  IF(ISNAN(rr)) CALL abort(__STAMP__,'HDG solver residual rr = NaN for CG iteration =', IntInfoOpt=iteration)
#if USE_MPI
  IF(MPIroot) converged=(rr.LT.AbortCrit2)

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

  CALL MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_PICLAS,iError)

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(3)  = MPIW8TimeField(3) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(3) = MPIW8CountField(3) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

#else
  converged=(rr.LT.AbortCrit2)
#endif /*USE_MPI*/
  IF(converged) THEN !converged
    TimeEndCG=PICLASTIME()
    CALL EvalResidual(1)
    CALL VectorDotProductRR(Norm_R2) !Z=V (function contains ALLREDUCE)
    IF(MPIroot) CALL DisplayConvergence(TimeEndCG-TimeStartCG, iteration, SQRT(Norm_R2))
    RETURN
  END IF !converged

  IF (MOD(iteration , MAX(INT(REAL(MaxIterCG)/REAL(OutIterCG)),1) ).EQ.0) THEN
    SWRITE(*,'(2(A,I0),2(A,G0))') 'CG solver reached ',iteration, ' of ',MaxIterCG, ' iterations with res = ',rr, ' > ',AbortCrit2
  END IF
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  IF(PrecondType.NE.0) THEN
    CALL ApplyPrecond(iVar,.FALSE.) ! ApplyPrecond(R,Z)
  ELSE
    DO SideID = 1, nSides
      HDG_Surf_N(SideID)%Z(iVar,:) = HDG_Surf_N(SideID)%R(iVar,:)
    END DO ! SideID = 1, nSides
  END IF
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/
  CALL VectorDotProductRZ(rz2)
#if USE_LOADBALANCE
  CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
  rz1 = rz2/rz1
  DO SideID = 1, nSides
    HDG_Surf_N(SideID)%V(iVar,:) = HDG_Surf_N(SideID)%Z(iVar,:) + rz1*HDG_Surf_N(SideID)%V(iVar,:)
  END DO ! SideID = 1, nSides
  rz1=rz2
#if USE_LOADBALANCE
  CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/
END DO ! iteration
SWRITE(*,*)'CG solver not converged in ',iteration, 'iterations!!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE CG_solver


!===================================================================================================================================
!> Set the global convergence properties of the HDG (CG) Solver and print then to StdOut)
!===================================================================================================================================
SUBROUTINE DisplayConvergence(ElapsedTime, iteration, Norm)
! MODULES
USE MOD_HDG_Vars      ,ONLY: HDGDisplayConvergence,HDGNorm,RunTime,RunTimePerIteration,iterationTotal,RunTimeTotal
USE MOD_Globals       ,ONLY: UNIT_StdOut
USE MOD_TimeDisc_Vars ,ONLY: iter,IterDisplayStep
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: ElapsedTime
INTEGER,INTENT(IN)  :: iteration
REAL,INTENT(IN)     :: Norm
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
RunTime = ElapsedTime
RunTimeTotal = RunTimeTotal + RunTime
IF(iteration.GT.0)THEN
  iterationTotal = iterationTotal + iteration
  RunTimePerIteration = RunTime/REAL(iteration)
ELSE
  RunTimePerIteration = 0.
END IF ! iteration.GT.0
HDGNorm = Norm

IF(HDGDisplayConvergence.AND.(MOD(iter,IterDisplayStep).EQ.0)) THEN
  WRITE(UNIT_StdOut,'(A,1X,I0,A,I0,A)')                '#iterations          :    ',iteration,' (',iterationTotal,' total)'
  WRITE(UNIT_StdOut,'(A,1X,ES25.14E3,A,ES25.14E3,A)')  'RunTime           [s]:',RunTime,' (',RunTimeTotal,' total)'
  WRITE(UNIT_StdOut,'(A,1X,ES25.14E3)')                'RunTime/iteration [s]:',RunTimePerIteration
  !WRITE(UNIT_StdOut,'(A,1X,ES16.7)')'RunTime/iteration/DOF[s]:',(TimeEndCG-TimeStartCG)/REAL(iteration*PP_nElems*nGP_vol)
  WRITE(UNIT_StdOut,'(A,1X,ES25.14E3)')                'Final Residual       :',HDGNorm
  WRITE(UNIT_StdOut,'(132("-"))')
END IF
END SUBROUTINE DisplayConvergence


!===================================================================================================================================
!> Computes Dot Product for vectors a and b: resu=a.b
!===================================================================================================================================
SUBROUTINE VectorDotProductRR(Resu)
! MODULES
USE MOD_Globals
USE MOD_PreProc
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars           ,ONLY: MPIW8TimeField,MPIW8CountField
#endif /*defined(MEASURE_MPI_WAIT)*/
USE MOD_HDG_Vars           ,ONLY: HDG_Surf_N
USE MOD_Mesh_Vars          ,ONLY: nSides
#if USE_MPI
USE MOD_Mesh_Vars          ,ONLY: nMPISides_YOUR
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!INTEGER,INTENT(IN):: dim1
!REAL,INTENT(IN)   :: A(dim1)
!REAL,INTENT(IN)   :: B(dim1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: SideID, endnSides
#if USE_MPI
REAL              :: ResuSend
#endif
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)   :: CounterStart,CounterEnd
REAL(KIND=8)      :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
#if USE_MPI
! not use MPI_YOUR sides for vector_dot_product!!!
endnSides=nSides-nMPIsides_YOUR
#else
endnSides=nSides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
Resu=0.
DO SideID = 1, endnSides
  Resu = Resu + DOT_PRODUCT(HDG_Surf_N(SideID)%R(1,:),HDG_Surf_N(SideID)%R(1,:))
END DO ! SideID = 1, nSides
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

#if USE_MPI
  ResuSend=Resu
  CALL MPI_ALLREDUCE(ResuSend,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(4)  = MPIW8TimeField(4) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(4) = MPIW8CountField(4) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

END SUBROUTINE VectorDotProductRR


!===================================================================================================================================
!> Performs matrix-vector multiplication for lambda system
!>   Parallel Mortar concept:
!>   1) MORTAR, BigToSmall: interpolate lambda from  big to small (small master sides)
!>   2) send lambda from master MPI sides to slave MPI sides (includes small mortar master sides)
!>   3) compute matrix-vector product locally on each proc, in mv array
!>   4) call mask_MPIsides: send  mv contribution from slave MPI sides to master MPI sides and add to master MPI sides
!>   5) MORTAR, SmallToBig: add contribution of finalized small mortar sides to big mortar, via Transpose of interpolation operator
!===================================================================================================================================
SUBROUTINE MatVec(iVar,DoVZ)
! MODULES
USE MOD_Globals
USE MOD_DG_Vars            ,ONLY: N_DG_Mapping
USE MOD_HDG_Vars           ,ONLY: nGP_face,nDirichletBCSides,DirichletBC,SetZeroPotentialDOF
USE MOD_HDG_Vars           ,ONLY: HDG_Surf_N,HDG_Vol_N
USE MOD_Mesh_Vars          ,ONLY: nSides, SideToElem, ElemToSide, nMPIsides_YOUR,N_SurfMesh, offSetElem
USE MOD_FillMortar_HDG     ,ONLY: BigToSmallMortar_HDG,SmallToBigMortar_HDG
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI                ,ONLY: StartReceiveMPISurfDataType,StartSendMPISurfDataType,FinishExchangeMPISurfDataType, Mask_MPIsides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBSplitTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
USE MOD_Interpolation_Vars ,ONLY: NMax,PREF_VDM
USE MOD_ChangeBasis        ,ONLY: ChangeBasis2D
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(INOUT) :: lambda(nGP_face, nSides)
LOGICAL,INTENT(IN),OPTIONAL :: DoVZ
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!REAL,INTENT(INOUT) :: mv(nGP_face, nSides)
INTEGER, INTENT(IN),OPTIONAL::iVar
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: firstSideID, lastSideID,NSideMin,Nloc,SideID2
INTEGER :: BCsideID,SideID, ElemID, locSideID
INTEGER :: jLocSide,jSideID(6)
#if USE_LOADBALANCE
REAL    :: tLBStart
#endif /*USE_LOADBALANCE*/
#if PP_nVar==1
INTEGER           :: dummy
#endif /*PP_nVar==1*/
REAL,DIMENSION(1:nGP_face(NMax)) :: lambdatmp,mvtmp
!===================================================================================================================================
mvtmp = 0.
lambdatmp = 0.

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
CALL BigToSmallMortar_HDG(1,DoVZ) ! lambda (DoVZ=F) or V (DoVZ=T)
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

#if USE_MPI
!CALL abort(__STAMP__,'not implemented')
!CALL StartReceiveMPIData(1,lambda,1,nSides, RecRequest_U,SendID=1) ! Receive YOUR
!CALL StartSendMPIData(   1,lambda,1,nSides,SendRequest_U,SendID=1) ! Send MINE
CALL StartReceiveMPISurfDataType(RecRequest_U , 1 , 2)
CALL StartSendMPISurfDataType(  SendRequest_U , 1 , MERGE(5,2,DoVZ), iVar) ! lambda (DoVZ=F) or V (DoVZ=T)
#endif /*USE_MPI*/


#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
firstSideID = 1
lastSideID = nSides-nMPIsides_YOUR

IF(DoVZ)THEN
  DO SideID=1,nSides
    HDG_Surf_N(SideID)%Z=0.
  END DO ! SideID=1,nSides
ELSE
  DO SideID=1,nSides
    HDG_Surf_N(SideID)%mv=0.
  END DO ! SideID=1,nSides
END IF ! DoVZ

DO SideID=firstSideID,lastSideID

  NSideMin = N_SurfMesh(SideID)%NSideMin
  !master element
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  IF(locSideID.NE.-1)THEN
    IF(DoVZ)THEN
      lambdatmp(1:nGP_face(NSideMin)) = HDG_Surf_N(SideID)%V(iVar,:)
    ELSE
      lambdatmp(1:nGP_face(NSideMin)) = HDG_Surf_N(SideID)%lambda(iVar,:)
    END IF ! DoVZ
    ElemID     = SideToElem(S2E_ELEM_ID,SideID)
    Nloc       = N_DG_Mapping(2,ElemID+offSetElem)
    IF(Nloc.GT.NSideMin)THEN
      ! From low to high
      CALL ChangeBasis2D(1, NSideMin, Nloc, PREF_VDM(NSideMin,Nloc)%Vdm , lambdatmp(1:nGP_face(NSideMin)), lambdatmp(1:nGP_face(Nloc)))
    END IF ! Nloc.GT.NSideMin
    jSideID(:) = ElemToSide(E2S_SIDE_ID,:,ElemID)
    DO jLocSide = 1,6
      SideID2 = jSideID(jLocSide)
      NSideMin = N_SurfMesh(SideID2)%NSideMin
      CALL DGEMV('N',nGP_face(Nloc),nGP_face(Nloc),1., &
                        HDG_Vol_N(ElemID)%Smat(:,:,jLocSide,locSideID), nGP_face(Nloc), &
                        lambdatmp(1:nGP_face(Nloc)),1,0.,& ! !: add to mv, 0: set mv
                        mvtmp(1:nGP_face(Nloc)),1)
      IF(Nloc.GT.NSideMin)THEN
        ! From high to low
        CALL ChangeBasis2D(1, Nloc, NSideMin, TRANSPOSE(PREF_VDM(NSideMin,Nloc)%Vdm) , mvtmp(1:nGP_face(Nloc)), mvtmp(1:nGP_face(NSideMin)))
      END IF ! Nloc.GT.NSideMin
      IF(DoVZ)THEN
        HDG_Surf_N(SideID2)%Z(iVar,:)  = HDG_Surf_N(SideID2)%Z(iVar,:)  + mvtmp(1:nGP_face(NSideMin))
      ELSE
        HDG_Surf_N(SideID2)%mv(iVar,:) = HDG_Surf_N(SideID2)%mv(iVar,:) + mvtmp(1:nGP_face(NSideMin))
      END IF ! DoVZ
    END DO !jLocSide
  END IF !locSideID.NE.-1

  NSideMin = N_SurfMesh(SideID)%NSideMin
  ! neighbour element
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  IF(locSideID.NE.-1)THEN
    IF(DoVZ)THEN
      lambdatmp(1:nGP_face(NSideMin)) = HDG_Surf_N(SideID)%V(iVar,:)
    ELSE
      lambdatmp(1:nGP_face(NSideMin)) = HDG_Surf_N(SideID)%lambda(iVar,:)
    END IF ! DoVZ
    ElemID     = SideToElem(S2E_NB_ELEM_ID,SideID)
    Nloc       = N_DG_Mapping(2,ElemID+offSetElem)
    jSideID(:) = ElemToSide(E2S_SIDE_ID,:,ElemID)
    IF(Nloc.GT.NSideMin)THEN
      ! From low to high
      CALL ChangeBasis2D(1, NSideMin, Nloc, PREF_VDM(NSideMin,Nloc)%Vdm , lambdatmp(1:nGP_face(NSideMin)), lambdatmp(1:nGP_face(Nloc)))
    END IF ! Nloc.GT.NSideMin
    DO jLocSide = 1,6
      SideID2 = jSideID(jLocSide)
      NSideMin = N_SurfMesh(SideID2)%NSideMin
      CALL DGEMV('N',nGP_face(Nloc),nGP_face(Nloc),1., &
                        HDG_Vol_N(ElemID)%Smat(:,:,jLocSide,locSideID), nGP_face(Nloc), &
                        lambdatmp(1:nGP_face(Nloc)),1,0.,& ! !: add to mv, 0: set mv
                        mvtmp(1:nGP_face(Nloc)),1)
      IF(Nloc.GT.NSideMin)THEN
        ! From high to low
        CALL ChangeBasis2D(1, Nloc, NSideMin, TRANSPOSE(PREF_VDM(NSideMin,Nloc)%Vdm) , mvtmp(1:nGP_face(Nloc)), mvtmp(1:nGP_face(NSideMin)))
      END IF ! Nloc.GT.NSideMin
      IF(DoVZ)THEN
        HDG_Surf_N(SideID2)%Z(iVar,:)  = HDG_Surf_N(SideID2)%Z(iVar,:)  + mvtmp(1:nGP_face(NSideMin))
      ELSE
        HDG_Surf_N(SideID2)%mv(iVar,:) = HDG_Surf_N(SideID2)%mv(iVar,:) + mvtmp(1:nGP_face(NSideMin))
      END IF ! DoVZ
    END DO !jLocSide
  END IF !locSideID.NE.-1
  !add mass matrix

END DO ! SideID=1,nSides
!SWRITE(*,*)'DEBUG---------------------------------------------------------'
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

#if USE_MPI
! Finish lambda communication
CALL FinishExchangeMPISurfDataType(SendRequest_U, RecRequest_U, 1, MERGE(5,2,DoVZ), iVar) ! lambda (DoVZ=F) or V (DoVZ=T)

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
firstSideID=nSides-nMPIsides_YOUR+1
lastSideID =nSides
DO SideID=firstSideID,lastSideID
  NSideMin = N_SurfMesh(SideID)%NSideMin
  !master element
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  IF(locSideID.NE.-1)THEN
    IF(DoVZ)THEN
      lambdatmp(1:nGP_face(NSideMin)) = HDG_Surf_N(SideID)%V(iVar,:)
    ELSE
      lambdatmp(1:nGP_face(NSideMin)) = HDG_Surf_N(SideID)%lambda(iVar,:)
    END IF ! DoVZ
    ElemID     = SideToElem(S2E_ELEM_ID,SideID)
    Nloc       = N_DG_Mapping(2,ElemID+offSetElem)
    IF(Nloc.GT.NSideMin)THEN
      ! From low to high
      CALL ChangeBasis2D(1, NSideMin, Nloc, PREF_VDM(NSideMin,Nloc)%Vdm , lambdatmp(1:nGP_face(NSideMin)), lambdatmp(1:nGP_face(Nloc)))
    END IF ! Nloc.GT.NSideMin
    jSideID(:) = ElemToSide(E2S_SIDE_ID,:,ElemID)
    DO jLocSide = 1,6
      SideID2 = jSideID(jLocSide)
      NSideMin = N_SurfMesh(SideID2)%NSideMin
      CALL DGEMV('N',nGP_face(Nloc),nGP_face(Nloc),1., &
                        HDG_Vol_N(ElemID)%Smat(:,:,jLocSide,locSideID), nGP_face(Nloc), &
                        lambdatmp(1:nGP_face(Nloc)),1,0.,& ! !: add to mv, 0: set mv
                        mvtmp(1:nGP_face(Nloc)),1)
      IF(Nloc.GT.NSideMin)THEN
        ! From high to low
        CALL ChangeBasis2D(1, Nloc, NSideMin, TRANSPOSE(PREF_VDM(NSideMin,Nloc)%Vdm) , mvtmp(1:nGP_face(Nloc)), mvtmp(1:nGP_face(NSideMin)))
      END IF ! Nloc.GT.NSideMin
      IF(DoVZ)THEN
        HDG_Surf_N(SideID2)%Z(iVar,:)  = HDG_Surf_N(SideID2)%Z(iVar,:)  + mvtmp(1:nGP_face(NSideMin))
      ELSE
        HDG_Surf_N(SideID2)%mv(iVar,:) = HDG_Surf_N(SideID2)%mv(iVar,:) + mvtmp(1:nGP_face(NSideMin))
      END IF ! DoVZ
    END DO !jLocSide
  END IF !locSideID.NE.-1

  NSideMin = N_SurfMesh(SideID)%NSideMin
  ! neighbour element
  locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  IF(locSideID.NE.-1)THEN
    IF(DoVZ)THEN
      lambdatmp(1:nGP_face(NSideMin)) = HDG_Surf_N(SideID)%V(iVar,:)
    ELSE
      lambdatmp(1:nGP_face(NSideMin)) = HDG_Surf_N(SideID)%lambda(iVar,:)
    END IF ! DoVZ
    ElemID     = SideToElem(S2E_NB_ELEM_ID,SideID)
    Nloc       = N_DG_Mapping(2,ElemID+offSetElem)
    jSideID(:) = ElemToSide(E2S_SIDE_ID,:,ElemID)
    IF(Nloc.GT.NSideMin)THEN
      ! From low to high
      CALL ChangeBasis2D(1, NSideMin, Nloc, PREF_VDM(NSideMin,Nloc)%Vdm , lambdatmp(1:nGP_face(NSideMin)), lambdatmp(1:nGP_face(Nloc)))
    END IF ! Nloc.GT.NSideMin
    DO jLocSide = 1,6
      SideID2 = jSideID(jLocSide)
      NSideMin = N_SurfMesh(SideID2)%NSideMin
      CALL DGEMV('N',nGP_face(Nloc),nGP_face(Nloc),1., &
                        HDG_Vol_N(ElemID)%Smat(:,:,jLocSide,locSideID), nGP_face(Nloc), &
                        lambdatmp(1:nGP_face(Nloc)),1,0.,& ! !: add to mv, 0: set mv
                        mvtmp(1:nGP_face(Nloc)),1)
      IF(Nloc.GT.NSideMin)THEN
        ! From high to low
        CALL ChangeBasis2D(1, Nloc, NSideMin, TRANSPOSE(PREF_VDM(NSideMin,Nloc)%Vdm) , mvtmp(1:nGP_face(Nloc)), mvtmp(1:nGP_face(NSideMin)))
      END IF ! Nloc.GT.NSideMin
      IF(DoVZ)THEN
        HDG_Surf_N(SideID2)%Z(iVar,:)  = HDG_Surf_N(SideID2)%Z(iVar,:)  + mvtmp(1:nGP_face(NSideMin))
      ELSE
        HDG_Surf_N(SideID2)%mv(iVar,:) = HDG_Surf_N(SideID2)%mv(iVar,:) + mvtmp(1:nGP_face(NSideMin))
      END IF ! DoVZ
    END DO !jLocSide
  END IF !locSideID.NE.-1
  !add mass matrix
END DO ! SideID=1,nSides
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/
IF (DoVZ) THEN
  CALL Mask_MPIsides('Z',iVar)
ELSE
  CALL Mask_MPIsides('mv',iVar)
END IF
#endif /*USE_MPI*/

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
CALL SmallToBigMortar_HDG(iVar,MERGE(1, 2, DoVZ)) ! CALL SmallToBigMortar_HDG(1,mv)

#if (PP_nVar!=1)
IF (iVar.EQ.4) THEN
#endif

  !set mv on Dirichlet BC to zero!
  IF(DoVZ)THEN
    DO BCsideID=1,nDirichletBCSides
      HDG_Surf_N(DirichletBC(BCsideID))%Z(iVar,:) =0.
    END DO ! SideID=1,nSides
  ELSE
    DO BCsideID=1,nDirichletBCSides
      HDG_Surf_N(DirichletBC(BCsideID))%mv(iVar,:)=0.
    END DO ! SideID=1,nSides
  END IF ! DoVZ

  ! Set potential to zero
  IF(SetZeroPotentialDOF)THEN
    IF(DoVZ)THEN
      HDG_Surf_N(1)%Z(iVar,1)  = 0.
    ELSe
      HDG_Surf_N(1)%mv(iVar,1) = 0.
    END IF ! DoVZ
  END IF ! SetZeroPotentialDOF

#if (PP_nVar!=1)
END IF
#endif

#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

! Suppress compiler warning
RETURN
#if PP_nVar==1
dummy=iVar
#endif /*PP_nVar==1*/

END SUBROUTINE MatVec


!===================================================================================================================================
!>
!===================================================================================================================================
SUBROUTINE EvalResidual(iVar)
! MODULES
USE MOD_Globals
USE MOD_HDG_Vars  ,ONLY: nDirichletBCSides,DirichletBC,SetZeroPotentialDOF
USE MOD_HDG_Vars  ,ONLY: HDG_Surf_N
USE MOD_Mesh_Vars ,ONLY: nSides
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL, INTENT(IN)    :: RHS(nGP_face,nSides)
!REAL, INTENT(INOUT) :: lambda(nGP_face,nSides)
INTEGER, INTENT(IN),OPTIONAL::iVar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!REAL, INTENT(OUT)   :: R(nGP_face,nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL                :: mv(nGP_face,nSides)
INTEGER             :: SideID,BCsideID
!===================================================================================================================================
! Set mv
IF(PRESENT(iVar)) THEN
  CALL MatVec(iVar,.FALSE.)
ELSE
  CALL MatVec(1,.FALSE.)
END IF

! R = RHS - mv
DO SideID = 1, nSides
  HDG_Surf_N(SideID)%R(iVar,:) = HDG_Surf_N(SideID)%RHS_face(iVar,:) - HDG_Surf_N(SideID)%mv(iVar,:)
END DO ! SideID = 1, nSides

!set mv on Dirichlet BC to zero!
#if (PP_nVar!=1)
IF (iVar.EQ.4) THEN
#endif
  ! Dirichlet BCs
  DO BCsideID=1,nDirichletBCSides
    HDG_Surf_N(DirichletBC(BCsideID))%R(iVar,:)=0.
  END DO ! SideID=1,nSides

  ! Set residual to zero
  IF(SetZeroPotentialDOF) HDG_Surf_N(1)%R(iVar,1) = 0.
#if (PP_nVar!=1)
END IF
#endif
END SUBROUTINE EvalResidual


!===================================================================================================================================
!> Apply the block-diagonal preconditioner for the lambda system
!===================================================================================================================================
SUBROUTINE ApplyPrecond(iVar,DoV)
! MODULES
USE MOD_Globals
USE MOD_HDG_Vars  ,ONLY: nGP_face,PrecondType,HDG_Surf_N,MaskedSide
USE MOD_Mesh_Vars ,ONLY: nSides, nMPIsides_YOUR,N_SurfMesh
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL,INTENT(IN) :: R(nGP_face, nSides)
LOGICAL,INTENT(IN)           :: DoV
INTEGER, INTENT(IN),OPTIONAL :: iVar
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!REAL,INTENT(INOUT) :: V(nGP_face, nSides)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: firstSideID, lastSideID, SideID, igf, NSideMin
!===================================================================================================================================
firstSideID = 1
lastSideID = nSides-nMPIsides_YOUR

IF(DoV)THEN

  ! Do V
  SELECT CASE(PrecondType)
  CASE(0)
    ! do nothing, should not be called
  CASE(1) !apply side-block SPD Preconditioner matrix, already Cholesky decomposed
    DO SideID=firstSideID,lastSideID
      IF(MaskedSide(sideID).NE.0) THEN
        HDG_Surf_N(SideID)%V(iVar,:)=0.
      ELSE
        NSideMin = N_SurfMesh(SideID)%NSideMin
        ! solve the preconditioner linear system
        CALL solveSPD(nGP_face(NSideMin),HDG_Surf_N(SideID)%Precond(:,:),1,HDG_Surf_N(SideID)%R(iVar,:), HDG_Surf_N(SideID)%V(iVar,:))
      END IF !maskedSide
    END DO ! SideID=1,nSides
  CASE(2)
    DO SideID=firstSideID,lastSideID
      IF(MaskedSide(sideID).NE.0) THEN
        HDG_Surf_N(SideID)%V(iVar,:)=0.
      ELSE
        NSideMin = N_SurfMesh(SideID)%NSideMin
        ! apply inverse of diagonal preconditioned
        DO igf = 1, nGP_face(NSideMin)
          HDG_Surf_N(SideID)%V(iVar,igf) = HDG_Surf_N(SideID)%InvPrecondDiag(igf)*HDG_Surf_N(SideID)%R(iVar,igf)
        END DO ! igf
      END IF !maskedSide
    END DO ! SideID=1,nSides
  END SELECT ! PrecondType

ELSE

  ! Do Z
  SELECT CASE(PrecondType)
  CASE(0)
    ! do nothing, should not be called
  CASE(1) !apply side-block SPD Preconditioner matrix, already Cholesky decomposed
    DO SideID=firstSideID,lastSideID
      IF(MaskedSide(sideID).NE.0) THEN
        HDG_Surf_N(SideID)%Z(iVar,:)=0.
      ELSE
        NSideMin = N_SurfMesh(SideID)%NSideMin
        ! solve the preconditioner linear system
        CALL solveSPD(nGP_face(NSideMin),HDG_Surf_N(SideID)%Precond(:,:),1,HDG_Surf_N(SideID)%R(iVar,:), HDG_Surf_N(SideID)%Z(iVar,:))
      END IF !maskedSide
    END DO ! SideID=1,nSides
  CASE(2)
    DO SideID=firstSideID,lastSideID
      IF(MaskedSide(sideID).NE.0) THEN
        HDG_Surf_N(SideID)%Z(iVar,:)=0.
      ELSE
        NSideMin = N_SurfMesh(SideID)%NSideMin
        ! apply inverse of diagonal preconditioned
        DO igf = 1, nGP_face(NSideMin)
          HDG_Surf_N(SideID)%Z(iVar,igf) = HDG_Surf_N(SideID)%InvPrecondDiag(igf)*HDG_Surf_N(SideID)%R(iVar,igf)
        END DO ! igf
      END IF !maskedSide
    END DO ! SideID=1,nSides
  END SELECT ! PrecondType

END IF ! DoV

END SUBROUTINE ApplyPrecond


!===================================================================================================================================
!> Solve a symmetrical positive definite linear system of dimension dims
!===================================================================================================================================
SUBROUTINE solveSPD(dimA,A,nRHS,RHS, X)
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN):: dimA
REAL,INTENT(IN)   :: A(dimA, dimA)
INTEGER,INTENT(IN):: nRHS
REAL,INTENT(IN)   :: RHS(dimA,nRHS)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT):: X(dimA,nRHS)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: lapack_info
!===================================================================================================================================
X = RHS
CALL DPOTRS('U',dimA,nRHS,A,dimA,X,dimA,lapack_info)
!IF (lapack_info .NE. 0) THEN
!  STOP 'LAPACK ERROR IN SOLVE CHOLESKY!'
!END IF
END SUBROUTINE solveSPD


!===================================================================================================================================
!> Computes Dot Product for vectors a and b: resu=a.b
!===================================================================================================================================
SUBROUTINE VectorDotProduct(dim1,A,B,Resu)
! MODULES
USE MOD_Globals
USE MOD_PreProc
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars           ,ONLY: MPIW8TimeField,MPIW8CountField
#endif /*defined(MEASURE_MPI_WAIT)*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN):: dim1
REAL,INTENT(IN)   :: A(dim1)
REAL,INTENT(IN)   :: B(dim1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
#if USE_MPI
REAL              :: ResuSend
#endif
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)   :: CounterStart,CounterEnd
REAL(KIND=8)      :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================

#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
Resu=0.
DO i=1,dim1
  Resu=Resu + A(i)*B(i)
END DO
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

#if USE_MPI
  ResuSend=Resu
  CALL MPI_ALLREDUCE(ResuSend,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(4)  = MPIW8TimeField(4) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(4) = MPIW8CountField(4) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

END SUBROUTINE VectorDotProduct


!===================================================================================================================================
!> Computes Dot Product for vectors a and b: resu=a.b
!===================================================================================================================================
SUBROUTINE VectorDotProductRZ(Resu)
! MODULES
USE MOD_Globals
USE MOD_PreProc
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars           ,ONLY: MPIW8TimeField,MPIW8CountField
#endif /*defined(MEASURE_MPI_WAIT)*/
USE MOD_HDG_Vars           ,ONLY: HDG_Surf_N
USE MOD_Mesh_Vars          ,ONLY: nSides
#if USE_MPI
USE MOD_Mesh_Vars          ,ONLY: nMPISides_YOUR
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!INTEGER,INTENT(IN):: dim1
!REAL,INTENT(IN)   :: A(dim1)
!REAL,INTENT(IN)   :: B(dim1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: SideID, endnSides
#if USE_MPI
REAL              :: ResuSend
#endif
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)   :: CounterStart,CounterEnd
REAL(KIND=8)      :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
#if USE_MPI
! not use MPI_YOUR sides for vector_dot_product!!!
endnSides=nSides-nMPIsides_YOUR
#else
endnSides=nSides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
Resu=0.
DO SideID = 1, endnSides
  Resu = Resu + DOT_PRODUCT(HDG_Surf_N(SideID)%R(1,:),HDG_Surf_N(SideID)%Z(1,:))
END DO ! SideID = 1, nSides
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

#if USE_MPI
  ResuSend=Resu
  CALL MPI_ALLREDUCE(ResuSend,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(4)  = MPIW8TimeField(4) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(4) = MPIW8CountField(4) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

END SUBROUTINE VectorDotProductRZ


!===================================================================================================================================
!> Computes Dot Product for vectors a and b: resu=a.b
!===================================================================================================================================
SUBROUTINE VectorDotProductRV(Resu)
! MODULES
USE MOD_Globals
USE MOD_PreProc
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars           ,ONLY: MPIW8TimeField,MPIW8CountField
#endif /*defined(MEASURE_MPI_WAIT)*/
USE MOD_HDG_Vars           ,ONLY: HDG_Surf_N
USE MOD_Mesh_Vars          ,ONLY: nSides
#if USE_MPI
USE MOD_Mesh_Vars          ,ONLY: nMPISides_YOUR
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!INTEGER,INTENT(IN):: dim1
!REAL,INTENT(IN)   :: A(dim1)
!REAL,INTENT(IN)   :: B(dim1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: SideID, endnSides
#if USE_MPI
REAL              :: ResuSend
#endif
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)   :: CounterStart,CounterEnd
REAL(KIND=8)      :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
#if USE_MPI
! not use MPI_YOUR sides for vector_dot_product!!!
endnSides=nSides-nMPIsides_YOUR
#else
endnSides=nSides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
Resu=0.
DO SideID = 1, endnSides
  Resu = Resu + DOT_PRODUCT(HDG_Surf_N(SideID)%R(1,:),HDG_Surf_N(SideID)%V(1,:))
END DO ! SideID = 1, nSides
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

#if USE_MPI
  ResuSend=Resu
  CALL MPI_ALLREDUCE(ResuSend,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(4)  = MPIW8TimeField(4) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(4) = MPIW8CountField(4) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

END SUBROUTINE VectorDotProductRV


!===================================================================================================================================
!> Computes Dot Product for vectors a and b: resu=a.b
!===================================================================================================================================
SUBROUTINE VectorDotProductVZ(Resu)
! MODULES
USE MOD_Globals
USE MOD_PreProc
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers ,ONLY: LBStartTime,LBPauseTime
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
USE MOD_MPI_Vars           ,ONLY: MPIW8TimeField,MPIW8CountField
#endif /*defined(MEASURE_MPI_WAIT)*/
USE MOD_HDG_Vars           ,ONLY: HDG_Surf_N
USE MOD_Mesh_Vars          ,ONLY: nSides
#if USE_MPI
USE MOD_Mesh_Vars          ,ONLY: nMPISides_YOUR
#endif
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!INTEGER,INTENT(IN):: dim1
!REAL,INTENT(IN)   :: A(dim1)
!REAL,INTENT(IN)   :: B(dim1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: Resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: SideID,endnSides
#if USE_MPI
REAL              :: ResuSend
#endif
#if USE_LOADBALANCE
REAL              :: tLBStart
#endif /*USE_LOADBALANCE*/
#if defined(MEASURE_MPI_WAIT)
INTEGER(KIND=8)   :: CounterStart,CounterEnd
REAL(KIND=8)      :: Rate
#endif /*defined(MEASURE_MPI_WAIT)*/
!===================================================================================================================================
#if USE_MPI
! not use MPI_YOUR sides for vector_dot_product!!!
endnSides=nSides-nMPIsides_YOUR
#else
endnSides=nSides
#endif /*USE_MPI*/
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/
Resu=0.
DO SideID = 1, endnSides
  Resu = Resu + DOT_PRODUCT(HDG_Surf_N(SideID)%V(1,:),HDG_Surf_N(SideID)%Z(1,:))
END DO ! SideID = 1, nSides
#if USE_LOADBALANCE
CALL LBPauseTime(LB_DG,tLBStart) ! Pause/Stop time measurement
#endif /*USE_LOADBALANCE*/

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterStart)
#endif /*defined(MEASURE_MPI_WAIT)*/

#if USE_MPI
  ResuSend=Resu
  CALL MPI_ALLREDUCE(ResuSend,Resu,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_PICLAS,iError)
#endif

#if defined(MEASURE_MPI_WAIT)
  CALL SYSTEM_CLOCK(count=CounterEnd, count_rate=Rate)
  MPIW8TimeField(4)  = MPIW8TimeField(4) + REAL(CounterEnd-CounterStart,8)/Rate
  MPIW8CountField(4) = MPIW8CountField(4) + 1_8
#endif /*defined(MEASURE_MPI_WAIT)*/

END SUBROUTINE VectorDotProductVZ



#endif /*USE_HDG*/


END MODULE MOD_HDG_Tools
