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

MODULE MOD_PICDepo_Tools
!===================================================================================================================================
! MOD PIC Depo
!===================================================================================================================================
IMPLICIT NONE
PRIVATE
!===================================================================================================================================
INTERFACE DepositParticleOnNodes
  MODULE PROCEDURE DepositParticleOnNodes
END INTERFACE

INTERFACE CalcCellLocNodeVolumes
  MODULE PROCEDURE CalcCellLocNodeVolumes
END INTERFACE

INTERFACE ReadTimeAverage
  MODULE PROCEDURE ReadTimeAverage
END INTERFACE

INTERFACE beta
  MODULE PROCEDURE beta
END INTERFACE

INTERFACE DeBoor
  MODULE PROCEDURE DeBoor
END INTERFACE

INTERFACE DeBoorRef
  MODULE PROCEDURE DeBoorRef
END INTERFACE

PUBLIC:: DepositParticleOnNodes,CalcCellLocNodeVolumes,ReadTimeAverage,beta,DeBoor,DeBoorRef
!===================================================================================================================================

CONTAINS


SUBROUTINE DepositParticleOnNodes(Charge,PartPos,ElemID) 
!----------------------------------------------------------------------------------------------------------------------------------!
! Deposit the charge of a single particle on the nodes corresponding to the deposition method 'cell_volweight_mean', where the
! charge density is stored in NodeSourceExt, which is added to NodeSource in the standard deposition procedure.
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PICDepo_Vars             ,ONLY: NodeSourceExt
USE MOD_Particle_Mesh_Vars       ,ONLY: GEO
#if ((USE_HDG) && (PP_nVar==1))
USE MOD_TimeDisc_Vars            ,ONLY: dt,tAnalyzeDiff,tEndDiff
#endif
USE MOD_Eval_xyz                 ,ONLY: GetPositionInRefElem
#if USE_LOADBALANCE
USE MOD_LoadBalance_Timers       ,ONLY: LBStartTime,LBElemPauseTime
#endif /*USE_LOADBALANCE*/
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
REAL,INTENT(IN)     :: Charge
REAL,INTENT(IN)     :: PartPos(1:3)
INTEGER,INTENT(IN)  :: ElemID
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                             :: alpha1, alpha2, alpha3, TempPartPos(1:3)
REAL                             :: TSource(1:4)
#if USE_LOADBALANCE
REAL                             :: tLBStart
#endif /*USE_LOADBALANCE*/
#if !((USE_HDG) && (PP_nVar==1))
INTEGER, PARAMETER               :: SourceDim=1
LOGICAL, PARAMETER               :: doCalculateCurrentDensity=.TRUE.
#else
LOGICAL                          :: doCalculateCurrentDensity
INTEGER                          :: SourceDim
#endif
!===================================================================================================================================
! Check whether charge and current density have to be computed or just the charge density
#if ((USE_HDG) && (PP_nVar==1))
IF(ALMOSTEQUAL(dt,tAnalyzeDiff).OR.ALMOSTEQUAL(dt,tEndDiff))THEN
doCalculateCurrentDensity=.TRUE.
SourceDim=1
ELSE ! do not calculate current density
doCalculateCurrentDensity=.FALSE.
SourceDim=4
END IF
#endif
#if USE_LOADBALANCE
CALL LBStartTime(tLBStart) ! Start time measurement
#endif /*USE_LOADBALANCE*/

CALL GetPositionInRefElem(PartPos,TempPartPos(1:3),ElemID,ForceMode=.TRUE.)
TSource(:) = 0.0
IF(doCalculateCurrentDensity)THEN
TSource(1:3) = PartPos*Charge
END IF
TSource(4) = Charge

alpha1=0.5*(TempPartPos(1)+1.0)
alpha2=0.5*(TempPartPos(2)+1.0)
alpha3=0.5*(TempPartPos(3)+1.0)
ASSOCIATE( NodeID     => GEO%ElemToNodeID(:,ElemID) )
NodeSourceExt(:,NodeID(1)) = NodeSourceExt(:,NodeID(1))+(TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*(1-alpha3))
NodeSourceExt(:,NodeID(2)) = NodeSourceExt(:,NodeID(2))+(TSource(SourceDim:4)*  (alpha1)*(1-alpha2)*(1-alpha3))
NodeSourceExt(:,NodeID(3)) = NodeSourceExt(:,NodeID(3))+(TSource(SourceDim:4)*  (alpha1)*  (alpha2)*(1-alpha3))
NodeSourceExt(:,NodeID(4)) = NodeSourceExt(:,NodeID(4))+(TSource(SourceDim:4)*(1-alpha1)*  (alpha2)*(1-alpha3))
NodeSourceExt(:,NodeID(5)) = NodeSourceExt(:,NodeID(5))+(TSource(SourceDim:4)*(1-alpha1)*(1-alpha2)*  (alpha3))
NodeSourceExt(:,NodeID(6)) = NodeSourceExt(:,NodeID(6))+(TSource(SourceDim:4)*  (alpha1)*(1-alpha2)*  (alpha3))
NodeSourceExt(:,NodeID(7)) = NodeSourceExt(:,NodeID(7))+(TSource(SourceDim:4)*  (alpha1)*  (alpha2)*  (alpha3))
NodeSourceExt(:,NodeID(8)) = NodeSourceExt(:,NodeID(8))+(TSource(SourceDim:4)*(1-alpha1)*  (alpha2)*  (alpha3))
END ASSOCIATE
#if USE_LOADBALANCE
CALL LBElemPauseTime(ElemID,tLBStart) ! Split time measurement (Pause/Stop and Start again) and add time to ElemID
#endif /*USE_LOADBALANCE*/

END SUBROUTINE DepositParticleOnNodes


SUBROUTINE CalcCellLocNodeVolumes()
!===================================================================================================================================
!> Initialize sub-cell volumes around nodes
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars          ,ONLY: sJ, nElems
USE MOD_Interpolation_Vars ,ONLY: wGP, xGP, wBary
USE MOD_ChangeBasis        ,ONLY: ChangeBasis3D
USE MOD_PreProc            ,ONLY: PP_N
USE MOD_Basis              ,ONLY: InitializeVandermonde
USE MOD_PICDepo_Vars       ,ONLY: CellLocNodes_Volumes
USE MOD_Particle_Mesh_Vars ,ONLY: GEO
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL    :: Vdm_loc(0:1,0:PP_N), wGP_loc, xGP_loc(0:1), DetJac(1,0:1,0:1,0:1)
REAL    :: DetLocal(1,0:PP_N,0:PP_N,0:PP_N)
INTEGER :: j, k, l, iElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CellLocNodes_Volumes = 0.0
IF (PP_N.NE.1) THEN
  xGP_loc(0) = -0.5
  xGP_loc(1) = 0.5
  wGP_loc = 1.
  CALL InitializeVandermonde(PP_N,1,wBary,xGP,xGP_loc, Vdm_loc)
END IF
DO iElem = 1, nElems
  IF (PP_N.EQ.1) THEN
    wGP_loc = wGP(0)
    DO j=0, PP_N; DO k=0, PP_N; DO l=0, PP_N
      DetJac(1,j,k,l)=1./sJ(j,k,l,iElem)
    END DO; END DO; END DO
  ELSE
    DO j=0, PP_N; DO k=0, PP_N; DO l=0, PP_N
      DetLocal(1,j,k,l)=1./sJ(j,k,l,iElem)
    END DO; END DO; END DO
    CALL ChangeBasis3D(1,PP_N, 1, Vdm_loc, DetLocal(:,:,:,:),DetJac(:,:,:,:))
  END IF
  ASSOCIATE( NodeVolume => CellLocNodes_Volumes(:),  &
             NodeID     => GEO%ElemToNodeID(:,iElem) )
    NodeVolume(NodeID(1)) = NodeVolume(NodeID(1)) + DetJac(1,0,0,0)
    NodeVolume(NodeID(2)) = NodeVolume(NodeID(2)) + DetJac(1,1,0,0)
    NodeVolume(NodeID(3)) = NodeVolume(NodeID(3)) + DetJac(1,1,1,0)
    NodeVolume(NodeID(4)) = NodeVolume(NodeID(4)) + DetJac(1,0,1,0)
    NodeVolume(NodeID(5)) = NodeVolume(NodeID(5)) + DetJac(1,0,0,1)
    NodeVolume(NodeID(6)) = NodeVolume(NodeID(6)) + DetJac(1,1,0,1)
    NodeVolume(NodeID(7)) = NodeVolume(NodeID(7)) + DetJac(1,1,1,1)
    NodeVolume(NodeID(8)) = NodeVolume(NodeID(8)) + DetJac(1,0,1,1)
  END ASSOCIATE
END DO

END SUBROUTINE CalcCellLocNodeVolumes


SUBROUTINE ReadTimeAverage(FileName)
!===================================================================================================================================
! Read in ChargeDensity and save to PartSource
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_IO_HDF5
USE MOD_HDF5_Input,              ONLY:ReadArray,ReadAttribute,File_ID,OpenDataFile,CloseDataFile,DatasetExists
USE MOD_Particle_Vars,           ONLY:nSpecies
USE MOD_PICDepo_Vars,            ONLY:PartSource
USE MOD_Mesh_Vars,               ONLY:OffsetElem,nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE         :: U(:,:,:,:,:)
INTEGER                  :: iSpec, iElem, kk, ll, mm
INTEGER                  :: Rank
INTEGER                  :: nVars, iVar, N_HDF5
INTEGER,ALLOCATABLE      :: PartSourceToVar(:)
INTEGER(HID_T)           :: Dset_ID,FileSpace
INTEGER(HSIZE_T), DIMENSION(7)          :: Dims,DimsMax
LOGICAL                  :: SolutionExists
CHARACTER(LEN=255),ALLOCATABLE          :: VarNames(:)
CHARACTER(LEN=10)        :: strhelp
#if USE_MPI
REAL                     :: StartT,EndT
#endif /*USE_MPI*/
!===================================================================================================================================

  SWRITE(UNIT_stdOut,*)'Using TimeAverage as constant PartSource(4) from file:',TRIM(FileName)
#if USE_MPI
  StartT=MPI_WTIME()
#endif

  IF(MPIRoot) THEN
    IF(.NOT.FILEEXISTS(FileName))  CALL abort(__STAMP__, &
          'TimeAverage-File "'//TRIM(FileName)//'" does not exist',999,999.)
  END IF
  CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_WORLD)

  ! get attributes
  CALL DatasetExists(File_ID,'DG_Solution',SolutionExists)
  IF(MPIRoot)THEN
    IF(.NOT.SolutionExists)  CALL abort(&
      __STAMP__&
      ,'DG_Solution in TimeAverage-File "'//TRIM(FileName)//'" does not exist!')
  END IF
  CALL H5DOPEN_F(File_ID, 'DG_Solution', Dset_ID, iError)
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
  IF(MPIRoot)THEN
    IF(INT(Dims(Rank),4).NE.nGlobalElems)  CALL abort(&
      __STAMP__&
      ,' MeshSize and Size of TimeAverage-File "'//TRIM(FileName)//'" does not match!')
  END IF
  nVars=INT(Dims(1),4)
  ALLOCATE(VarNames(nVars))
  CALL ReadAttribute(File_ID,'VarNames',nVars,StrArray=VarNames)

  ALLOCATE(PartSourceToVar(nSpecies))
  PartSourceToVar=0
  DO iSpec=1,nSpecies
    WRITE(strhelp,'(I2.2)') iSpec
    DO iVar=1,nVars
      IF (VarNames(iVar).EQ.TRIM('ChargeDensity-Spec')//TRIM(strhelp)) THEN
        PartSourceToVar(iSpec)=iVar
        EXIT
      END IF
    END DO
  END DO
  IF (.NOT.ANY(PartSourceToVar.NE.0)) CALL abort(__STAMP__, &
    'No PartSource found in TimeAverage-File "'//TRIM(FileName)//'"!!!',999,999.)
  DEALLOCATE(VarNames)

  !-- read state
  ALLOCATE(U(nVars,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
  CALL ReadAttribute(File_ID,'N',1,IntegerScalar=N_HDF5)
  IF(N_HDF5.EQ.PP_N)THEN! No interpolation needed, read solution directly from file
    ! Associate construct for integer KIND=8 possibility
    ASSOCIATE (&
          nVars       => INT(nVars,IK)     ,&
          PP_N        => INT(PP_N,IK)      ,&
          PP_nElems   => INT(PP_nElems,IK) ,&
          OffsetElem  => INT(OffsetElem,IK) )
          CALL ReadArray('DG_Solution',5,(/nVars,PP_N+1_IK,PP_N+1_IK,PP_N+1_IK,PP_nElems/),OffsetElem,5,RealArray=U)
    END ASSOCIATE
  ELSE
    CALL abort(__STAMP__, &
          'N_HDF5.NE.PP_N !',999,999.)
  END IF
  CALL CloseDataFile()

  !-- save to PartSource
  PartSource(4,:,:,:,:)=0.
  DO iSpec=1,nSpecies
    IF (PartSourceToVar(iSpec).NE.0) THEN
      DO iElem=1,PP_nElems
        DO kk = 0, PP_N
          DO ll = 0, PP_N
            DO mm = 0, PP_N
              PartSource(4,mm,ll,kk,iElem)=PartSource(4,mm,ll,kk,iElem)+U(PartSourceToVar(iSpec),mm,ll,kk,iElem)
            END DO
          END DO
        END DO
      END DO
    END IF
  END DO
  DEALLOCATE(U)
  DEALLOCATE(PartSourceToVar)

#if USE_MPI
  EndT=MPI_WTIME()
  SWRITE(UNIT_stdOut,'(A,F0.3,A)',ADVANCE='YES')' Readin took  [',EndT-StartT,'s].'
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' DONE!'
#else
  SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')' DONE!'
#endif

END SUBROUTINE ReadTimeAverage


SUBROUTINE DeBoor(PosInd, aux, coord, results, dir)
!============================================================================================================================
! recursive function for evaluating a b-spline basis function
!============================================================================================================================
! use MODULES
USE MOD_PICDepo_Vars
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                          :: PosInd, dir
REAL, INTENT(IN)                             :: coord
REAL, INTENT(INOUT)                          :: aux(0:3), results
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,k,jL,jR
REAL                              :: hlp1,hlp2
!-----------------------------------------------------------------------------------------------------------------------------------
DO i = 0, 2
   DO k = 0, 2-i
      jL = PosInd - k
      jR = jL + 3 - i

      hlp1 = jR * BGMdeltas(dir) - coord
      hlp2 = coord - jL * BGMdeltas(dir)

      aux(k) = (hlp1 * aux(k+1) + hlp2 * aux(k)) / (hlp1+hlp2)
   ENDDO
ENDDO
results = aux(0)

END SUBROUTINE DeBoor


SUBROUTINE DeBoorRef(N_in,NKnots,Knots,Xi_in,UBspline)
!===================================================================================================================================
! DeBoor algorithms for uniform B-Splines
! rule of DeBoor algorithm
! N_i^0 (x) = 1 if x in [u_i, u_i+1); else 0
! N_i^n (x) = (x-u_i) /(u_i+n-u_i) N_i^n-1(x) + (u_i+n+1 - x)/(u_i+n+1 - u_i+1) N_i+1^n-1(x)
! this algorithm evaluates the complete 1D basis fuction, because certain knots can be reused
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
INTEGER,INTENT(IN)      :: N_in,Nknots
REAL,INTENT(IN)         :: Xi_in
REAL,INTENT(OUT)        :: UBspline(0:N_in)
REAL,INTENT(IN)         :: knots(0:Nknots) ! range of parameter
!REAL,INTENT(IN)         :: DXi is 2 for [-1,1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,n, last!,first
REAL                    :: tmpArray(0:NKnots-1,0:NKnots-1)
REAL                    :: sDxiN!,DxiN1
!===================================================================================================================================


! init, first layer (constant)
! zero order
n=0
tmpArray=0.
last=nknots-1
DO i=0,last
  IF((knots(i).LE.Xi_in).AND.(Xi_in.LT.knots(i+1)))THEN
    tmpArray(i,n)=1.0
  END IF ! select
END DO ! i

DO n=1,N_in
  last=last-1
  DO i=0,last
    ! standard
!    tmpArray(i,n) = (Xi_in-knots(i))/(knots(i+n)-knots(i))*tmpArray(i,n-1) &
!                  + (knots(i+n+1)-Xi_in)/(knots(i+n+1)-knots(i+1))*tmpArray(i+1,n-1)
    ! optimized
    sDxiN=0.5/REAL(n)
    tmpArray(i,n) = sDxiN*( (Xi_in-knots(i) )*tmparray(i,n-1)+(knots(i+n+1)-Xi_in)*tmpArray(i+1,n-1))
  END DO ! i
END DO ! n

! move back to correct range
UBSpline(0:N_in)=tmpArray(0:N_in,N_in)

END SUBROUTINE DeBoorRef


FUNCTION beta(z,w)
!============================================================================================================================
! calculates the beta function
!============================================================================================================================
! use MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
   IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
   REAL, INTENT(IN)                             :: w, z
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
   REAL                                          :: beta
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
   beta = GAMMA(z)*GAMMA(w)/GAMMA(z+w)
END FUNCTION beta


END MODULE MOD_PICDepo_Tools
