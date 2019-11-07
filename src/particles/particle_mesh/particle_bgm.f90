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

MODULE MOD_Particle_BGM
!===================================================================================================================================
!> Contains
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

PUBLIC::DefineParametersParticleBGM
PUBLIC::BuildBGM

CONTAINS

!==================================================================================================================================
!> Define parameters for particle backgroundmesh
!==================================================================================================================================
SUBROUTINE DefineParametersParticleBGM()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection('BGM')

! Background mesh init variables
CALL prms%CreateRealArrayOption('Part-FIBGMdeltas'&
  , 'Define the deltas for the cartesian Fast-Init-Background-Mesh.'//&
  ' They should be of the similar size as the smallest cells of the used mesh for simulation.'&
  , '1. , 1. , 1.')
CALL prms%CreateRealArrayOption('Part-FactorFIBGM'&
  , 'Factor with which the background mesh will be scaled.'&
  , '1. , 1. , 1.')

END SUBROUTINE DefineParametersParticleBGM


SUBROUTINE BuildBGM()
!===================================================================================================================================
!> computes the BGM-indices of an element and maps the number of element and which element to each BGM cell
!> BGM is only saved for compute-node-mesh + halo-region on shared memory
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars            ,ONLY: NodeCoords, nElems
USE MOD_Partilce_Periodic_BC ,ONLY: InitPeriodicBC
USE MOD_Particle_Mesh_Vars   ,ONLY: GEO
USE MOD_Particle_MPI_Vars    ,ONLY: SafetyFactor,halo_eps_velo,halo_eps,halo_eps2
USE MOD_Equation_Vars        ,ONLY: c
USE MOD_Particle_Vars        ,ONLY: manualtimestep
USE MOD_ReadInTools          ,ONLY: GetRealArray
#if !(USE_HDG)
USE MOD_CalcTimeStep         ,ONLY: CalcTimeStep
#endif /*USE_HDG*/
#if USE_MPI
USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared           ,ONLY: Allocate_Shared
USE MOD_PICDepo_Vars         ,ONLY: DepositionType, r_sf
#endif /*USE_MPI*/
USE MOD_ReadInTools          ,ONLY: PrintOption
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem
INTEGER,ALLOCATABLE            :: ElemToBGM(:,:)
REAL                           :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER                        :: iBGM, jBGM, kBGM
INTEGER                        :: BGMimax, BGMimin, BGMjmax, BGMjmin, BGMkmax, BGMkmin
INTEGER                        :: BGMCellXmax, BGMCellXmin, BGMCellYmax, BGMCellYmin, BGMCellZmax, BGMCellZmin
REAL                           :: deltaT
REAL                           :: globalDiag
#if USE_MPI
INTEGER,ALLOCATABLE            :: sendbuf(:,:,:),recvbuf(:,:,:)
INTEGER,ALLOCATABLE            :: offsetElemsInBGMCell(:,:,:)
INTEGER(KIND=MPI_ADDRESS_KIND) :: MPISharedSize
#endif
!===================================================================================================================================

!! Read parameter for FastInitBackgroundMesh (FIBGM)
GEO%FIBGMdeltas(1:3) = GETREALARRAY('Part-FIBGMdeltas',3,'1. , 1. , 1.')
GEO%FactorFIBGM(1:3) = GETREALARRAY('Part-FactorFIBGM',3,'1. , 1. , 1.')
GEO%FIBGMdeltas(1:3) = 1./GEO%FactorFIBGM(1:3) * GEO%FIBGMdeltas(1:3)

MPISharedSize = INT(6*nTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/6,nTotalElems/),ElemToBGM_Shared_Win,ElemToBGM_Shared)
CALL Allocate_Shared(MPISharedSize,(/6,nTotalElems/),BoundsOfElem_Shared_Win,BoundsOfElem_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToBGM_Shared_Win,IERROR)
CALL MPI_WIN_LOCK_ALL(0,BoundsOfElem_Shared_Win,IERROR)

firstElem=INT(REAL(myRank_Shared*nTotalElems)/REAL(nProcessors_Shared))+1
lastElem=INT(REAL((myRank_Shared+1)*nTotalElems)/REAL(nProcessors_Shared))
DO iElem = firstElem, lastElem
  firstNodeID=ElemInfo_Shared(ELEM_FirstNodeInd,iElem)
  nNodeIDs=ElemInfo_Shared(ELEM_LastNodeInd,iElem)-ElemInfo_Shared(ELEM_FirstNodeInd,iElem)

  xmin=MINVAL(NodeCoords_Shared(1,firstNodeID:firstNodeID+nNodeIDs))
  xmax=MAXVAL(NodeCoords_Shared(1,firstNodeID:firstNodeID+nNodeIDs))
  ymin=MINVAL(NodeCoords_Shared(2,firstNodeID:firstNodeID+nNodeIDs))
  ymax=MAXVAL(NodeCoords_Shared(2,firstNodeID:firstNodeID+nNodeIDs))
  zmin=MINVAL(NodeCoords_Shared(3,firstNodeID:firstNodeID+nNodeIDs))
  zmax=MAXVAL(NodeCoords_Shared(3,firstNodeID:firstNodeID+nNodeIDs))

  BoundsOfElem_Shared(1,iElem) = xmin
  BoundsOfElem_Shared(2,iElem) = xmax
  BoundsOfElem_Shared(3,iElem) = ymin
  BoundsOfElem_Shared(4,iElem) = ymax
  BoundsOfElem_Shared(5,iElem) = zmin
  BoundsOfElem_Shared(6,iElem) = zmax

  ElemToBGM_Shared(1,iElem) = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
  ElemToBGM_Shared(2,iElem) = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
  ElemToBGM_Shared(3,iElem) = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
  ElemToBGM_Shared(4,iElem) = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
  ElemToBGM_Shared(5,iElem) = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
  ElemToBGM_Shared(6,iElem) = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
END DO ! iElem = 1, nElems
CALL MPI_WIN_SYNC(ElemToBGM_Shared_Win,IERROR)
CALL MPI_WIN_SYNC(BoundsOfElem_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,IERROR)

!CALL InitPeriodicBC()

! deallocate stuff // required for dynamic load balance
#if USE_MPI
IF (ALLOCATED(GEO%FIBGM)) THEN
  DO iBGM=GEO%FIBGMimin,GEO%FIBGMimax
    DO jBGM=GEO%FIBGMjmin,GEO%FIBGMjmax
      DO kBGM=GEO%FIBGMkmin,GEO%FIBGMkmax
        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element)
        !SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%ShapeProcs)
        !SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs)
        !SDEALLOCATE(GEO%FIBGM(i,k,l)%SharedProcs)
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
  DEALLOCATE(GEO%FIBGM)
END IF
#endif /*USE_MPI*/

!--- Read Manual Time Step
useManualTimeStep = .FALSE.
ManualTimeStep = GETREAL('Particles-ManualTimeStep', '0.0')
IF (ManualTimeStep.GT.0.0) THEN
  useManualTimeStep=.True.
END IF
SafetyFactor  =GETREAL('Part-SafetyFactor','1.0')
halo_eps_velo =GETREAL('Particles-HaloEpsVelo','0')

IF (ManualTimeStep.EQ.0.0) THEN
#if !(USE_HDG)
  deltaT=CALCTIMESTEP()
#else
   CALL abort(&
__STAMP__&
, 'ManualTimeStep.EQ.0.0 -> ManualTimeStep is not defined correctly! Particles-ManualTimeStep = ',RealInfoOpt=ManualTimeStep)
#endif /*USE_HDG*/
ELSE
  deltaT=ManualTimeStep
END IF
IF (halo_eps_velo.EQ.0) halo_eps_velo = c
#if (PP_TimeDiscMethod==4 || PP_TimeDiscMethod==200 || PP_TimeDiscMethod==42 || PP_TimeDiscMethod==43)
IF (halo_eps_velo.EQ.c) THEN
   CALL abort(&
__STAMP__&
, 'halo_eps_velo.EQ.c -> Halo Eps Velocity for MPI not defined')
END IF
#endif
#if (PP_TimeDiscMethod==501) || (PP_TimeDiscMethod==502) || (PP_TimeDiscMethod==506)
halo_eps = RK_c(2)
DO iStage=2,nRKStages-1
  halo_eps = MAX(halo_eps,RK_c(iStage+1)-RK_c(iStage))
END DO
halo_eps = MAX(halo_eps,1.-RK_c(nRKStages))
CALL PrintOption('max. RKdtFrac','CALCUL.',RealOpt=halo_eps)
halo_eps = halo_eps*halo_eps_velo*deltaT*SafetyFactor !dt multiplied with maximum RKdtFrac
#else
halo_eps = halo_eps_velo*deltaT*SafetyFactor ! for RK too large
#endif

#if USE_MPI
! Check whether halo_eps is smaller than shape function radius
! e.g. 'shape_function', 'shape_function_1d', 'shape_function_cylindrical', 'shape_function_spherical', 'shape_function_simple'
IF(TRIM(DepositionType(1:MIN(14,LEN(TRIM(ADJUSTL(DepositionType)))))).EQ.'shape_function')THEN
  IF(halo_eps.LT.r_sf)THEN
    SWRITE(UNIT_stdOut,'(A)') ' halo_eps is smaller than shape function radius. Setting halo_eps=r_sf'
    halo_eps = halo_eps + r_sf
    CALL PrintOption('max. RKdtFrac','CALCUL.',RealOpt=halo_eps)
  END IF
END IF
#endif /*USE_MPI*/

! limit halo_eps to diagonal of bounding box
globalDiag = SQRT( (GEO%xmaxglob-GEO%xminglob)**2 &
                 + (GEO%ymaxglob-GEO%yminglob)**2 &
                 + (GEO%zmaxglob-GEO%zminglob)**2 )
IF(halo_eps.GT.globalDiag)THEN
  CALL PrintOption('unlimited halo distance','CALCUL.',RealOpt=halo_eps)
  SWRITE(UNIT_stdOut,'(A38)') ' |   limitation of halo distance  |    '
  halo_eps=globalDiag
END IF

halo_eps2=halo_eps*halo_eps
CALL PrintOption('halo distance','CALCUL.',RealOpt=halo_eps)

! initialize BGM min/max indeces using GEO min/max distances
#if USE_MPI
! enlarge BGM with halo region (all element outside of this region will be cut off)
BGMimax = INT((MIN(GEO%xmax_Shared+halo_eps,GEO%xmaxglob)-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
BGMimin = INT((MAX(GEO%xmin_Shared-halo_eps,GEO%xminglob)-GEO%xminglob)/GEO%FIBGMdeltas(1))-1
BGMjmax = INT((MIN(GEO%ymax_Shared+halo_eps,GEO%ymaxglob)-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
BGMjmin = INT((MAX(GEO%ymin_Shared-halo_eps,GEO%yminglob)-GEO%yminglob)/GEO%FIBGMdeltas(2))-1
BGMkmax = INT((MIN(GEO%zmax_Shared+halo_eps,GEO%zmaxglob)-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
BGMkmin = INT((MAX(GEO%zmin_Shared-halo_eps,GEO%zminglob)-GEO%zminglob)/GEO%FIBGMdeltas(3))-1
#else
BGMimax = INT((GEO%xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))+1
BGMimin = INT((GEO%xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))-1
BGMjmax = INT((GEO%ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))+1
BGMjmin = INT((GEO%ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))-1
BGMkmax = INT((GEO%zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))+1
BGMkmin = INT((GEO%zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))-1
#endif /*USE_MPI*/

! write function-local BGM indeces into global variables
GEO%FIBGMimax_Shared=BGMimax
GEO%FIBGMimin_Shared=BGMimin
GEO%FIBGMjmax_Shared=BGMjmax
GEO%FIBGMjmin_Shared=BGMjmin
GEO%FIBGMkmax_Shared=BGMkmax
GEO%FIBGMkmin_Shared=BGMkmin

ALLOCATE(GEO%FIBGM(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax))

! null number of element per BGM cell
DO kBGM = BGMkmin,BGMkmax
  DO jBGM = BGMjmin,BGMjmax
    DO iBGM = BGMimin,BGMimax
      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

!--- compute number of elements in each background cell
! allocated shared memory for nElems per BGM cell

! check which element is inside of compute-node domain (1), check which element is inside of compute-node halo (2)
! and which element is outside of compute-node domain (0)
! first do coarse check with BGM
ElemInfo_Shared(ElemInfoSize+1,firstElem:lastElem)=0
DO iElem = firstElem, lastElem
  BGMCellXmin = ElemToBGM_Shared(1,iElem)
  BGMCellXmax = ElemToBGM_Shared(2,iElem)
  BGMCellYmin = ElemToBGM_Shared(3,iElem)
  BGMCellYmax = ElemToBGM_Shared(4,iElem)
  BGMCellZmin = ElemToBGM_Shared(5,iElem)
  BGMCellZmax = ElemToBGM_Shared(6,iElem)
  ! add current element to number of BGM-elems
  DO iBGM = BGMCellXmin,BGMCellXmax
    IF(iBGM.LT.BGMimin) CYCLE
    IF(iBGM.GT.BGMimax) CYCLE
    DO jBGM = BGMCellYmin,BGMCellYmax
      IF(jBGM.LT.BGMjmin) CYCLE
      IF(jBGM.GT.BGMjmax) CYCLE
      DO kBGM = BGMCellZmin,BGMCellZmax
        IF(kBGM.LT.BGMkmin) CYCLE
        IF(kBGM.GT.BGMkmax) CYCLE
        !GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
        IF(iElem.GT.offsetElem_Shared+1 .AND. iElem.LE.offsetElem_Shared+nElems_Shared) THEN
          ElemInfo_Shared(ElemInfoSize+1,iElem)=1 ! compute-node element
        ELSE
          ElemInfo_Shared(ElemInfoSize+1,iElem)=2 ! halo element
        END IF
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem
CALL MPI_WIN_SYNC(ElemInfo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

! sum up potential halo elements and create correct offset mapping
nHaloElems = 0
ALLOCATE(offsetHaloElem(nTotalElems))
DO iElem = 1, nTotalElems
  IF (ElemInfo_Shared(ElemInfoSize+1,iElem).EQ.2) THEN
    nHaloElems = nHaloElems + 1
    offsetHaloElem(nHaloElems) = iElem
  END IF
END DO

! sum all MPI-side of compute-node and create correct offset mapping
nMPISides_Shared = 0
ALLOCATE(offsetMPISide_Shared(nTotalSides))
DO iSide = 1, nTotalSides
  IF (SideInfo_Shared(SideInfoSize+1,iSide).EQ.2) THEN
    nMPISides_Shared = nMPISides_Shared + 1
    offsetMPISide_Shared(nMPISides_Shared) = iSide
  END IF
END DO

! Distribute nHaloElements evenly on compute-node procs
IF (nHaloElems.GT.nProcessors_Shared)
  firstHaloElem=INT(REAL(myRank_Shared*nHaloElems)/REAL(nProcessors_Shared))+1
  lastHaloElem=INT(REAL((myRank_Shared+1)*nHaloElems)/REAL(nProcessors_Shared))
ELSE
  firstHaloElem = myRank_Shared + 1
  IF (myRank_Shared.LT.nHaloElems) THEN
    lastHaloElem = myRank_Shared + 1
  ELSE
    lastHaloElem = 0
  END IF
END IF

! do refined check:
! check the bounding box of each element in compute-node halo domain 
! against the bounding boxes of compute-node MPI-border faces' element
DO iElem = firstHaloElem, lastHaloElem
  ElemID = offsetHaloElem(iElem)
  ElemInsideHalo = .FALSE.
  DO iSide = 1, nMPISides_Shared
    SideID = offsetMPISide_Shared(iSide)
    IF ((BoundsOfElem_Shared(1,SideInfo_Shared(SIDE_ElemID,SideID))-halo_eps).GT.BoundsOfElem(2,ElemID)) CYCLE
    IF ((BoundsOfElem_Shared(2,SideInfo_Shared(SIDE_ElemID,SideID))+halo_eps).LT.BoundsOfElem(1,ElemID)) CYCLE
    IF ((BoundsOfElem_Shared(3,SideInfo_Shared(SIDE_ElemID,SideID))-halo_eps).GT.BoundsOfElem(4,ElemID)) CYCLE
    IF ((BoundsOfElem_Shared(4,SideInfo_Shared(SIDE_ElemID,SideID))+halo_eps).LT.BoundsOfElem(3,ElemID)) CYCLE
    IF ((BoundsOfElem_Shared(5,SideInfo_Shared(SIDE_ElemID,SideID))-halo_eps).GT.BoundsOfElem(6,ElemID)) CYCLE
    IF ((BoundsOfElem_Shared(6,SideInfo_Shared(SIDE_ElemID,SideID))+halo_eps).LT.BoundsOfElem(5,ElemID)) CYCLE
    ElemInsideHalo = .TRUE.
    EXIT
  END DO ! iSide = 1, nMPISides_Shared
  IF (.NOT.ElemInsideHalo) THEN
    ElemInfo_Shared(ElemInfoSize+1,ElemID)=0
  ELSE
    BGMCellXmin = ElemToBGM_Shared(1,ElemID)
    BGMCellXmax = ElemToBGM_Shared(2,ElemID)
    BGMCellYmin = ElemToBGM_Shared(3,ElemID)
    BGMCellYmax = ElemToBGM_Shared(4,ElemID)
    BGMCellZmin = ElemToBGM_Shared(5,ElemID)
    BGMCellZmax = ElemToBGM_Shared(6,ElemID)
    ! add current element to number of BGM-elems
    DO iBGM = BGMCellXmin,BGMCellXmax
      DO jBGM = BGMCellYmin,BGMCellYmax
        DO kBGM = BGMCellZmin,BGMCellZmax
          GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
        END DO ! kBGM
      END DO ! jBGM
    END DO ! iBGM
  END IF
END DO ! iElem = firstHaloElem, lastHaloElem
CALL MPI_WIN_SYNC(ElemInfo_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

!--- compute number of elements in each background cell
DO iElem = offsetElem+1, offsetElem+nElems
  BGMCellXmin = ElemToBGM_Shared(1,iElem)
  BGMCellXmax = ElemToBGM_Shared(2,iElem)
  BGMCellYmin = ElemToBGM_Shared(3,iElem)
  BGMCellYmax = ElemToBGM_Shared(4,iElem)
  BGMCellZmin = ElemToBGM_Shared(5,iElem)
  BGMCellZmax = ElemToBGM_Shared(6,iElem)
  ! add current element to number of BGM-elems
  DO iBGM = BGMCellXmin,BGMCellXmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO kBGM = BGMCellZmin,BGMCellZmax
        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem

! alternative nElem count with cycles
!DO iElem = firstElem, lastElem
!  IF (ElemInfo_Shared(ElemInfoSize+1,iElem).EQ.0) CYCLE
!  BGMCellXmin = ElemToBGM_Shared(1,iElem)
!  BGMCellXmax = ElemToBGM_Shared(2,iElem)
!  BGMCellYmin = ElemToBGM_Shared(3,iElem)
!  BGMCellYmax = ElemToBGM_Shared(4,iElem)
!  BGMCellZmin = ElemToBGM_Shared(5,iElem)
!  BGMCellZmax = ElemToBGM_Shared(6,iElem)
!  ! add current element to number of BGM-elems
!  DO iBGM = BGMCellXmin,BGMCellXmax
!    DO jBGM = BGMCellYmin,BGMCellYmax
!      DO kBGM = BGMCellZmin,BGMCellZmax
!        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
!      END DO ! kBGM
!    END DO ! jBGM
!  END DO ! iBGM
!END DO ! iElem

ALLOCATE(sendbuf(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax))
ALLOCATE(recvbuf(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax))
! find max nelems and offset in each BGM cell
DO iBGM = BGMimin,BGMimax
  DO jBGM = BGMjmin,BGMjmax
    DO kBGM = BGMkmin,BGMkmax
      sendbuf(iBGM,jBGM,kBGM)=GEO%FIBGM(iBGM,jBGM,kBGM)%nElem
      recvbuf(iBGM,jBGM,kBGM)=0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

ALLOCATE(offsetElemsInBGMCell(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax))
CALL MPI_EXSCAN(sendbuf(:,:,:),recvbuf(:,:,:),((BGMimax-BGMimin)+1)*((BGMjmax-BGMjmin)+1)*((BGMkmax-BGMkmin)+1) &
                ,MPI_INTEGER,MPI_SUM,MPI_COMM_SHARED,iError)
offsetElemsInBGMCell=recvbuf
DEALLOCATE(recvbuf)

! last proc of compute-node calculates total number of elements in each BGM-cell 
! after this loop sendbuf of last proc contains nElems per BGM cell
IF(myRank_Shared.EQ.nProcessors_Shared-1)THEN
  DO iBGM = BGMimin,BGMimax
    DO jBGM = BGMjmin,BGMjmax
      DO kBGM = BGMkmin,BGMkmax
        sendbuf(iBGM,jBGM,kBGM)=offsetElemsInBGMCell(iBGM,jBGM,kBGM)+GEO%FIBGM(iBGM,jBGM,kBGM)%nElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END IF

! allocated shared memory for nElems per BGM cell
MPISharedSize = INT(((BGMimax-BGMimin)+1)*((BGMjmax-BGMjmin)+1)*((BGMkmax-BGMkmin)+1),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/BGMimax-BGMimin+1,BGMjmax-BGMjmin+1,BGMkmax-BGMkmin+1/) &
                    ,FIBGM_nElem_Shared_Win,FIBGM_nElem_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGM_nElem_Shared_Win,IERROR)

! last proc of compute-node writes into shared memory to make nElems per BGM accessible for every proc
IF(myRank_Shared.EQ.nProcessors_Shared-1)THEN
  DO iBGM = BGMimin,BGMimax
    DO jBGM = BGMjmin,BGMjmax
      DO kBGM = BGMkmin,BGMkmax
        FIBGM_nElem_Shared(iBGM,jBGM,kBGM) = sendbuf(iBGM,jBGM,kBGM)
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END IF
DEALLOCATE(sendbuf)
CALL MPI_WIN_SYNC(FIBGM_nElem_Shared_Win,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

! allocate 1D array for mapping of BGM cell to Element indeces
MPISharedSize = INT((FIBGM_offsetElem_Shared(BGMimax,BGMjmax,BGMkmax)+FIBGM_nElem_Shared(BGMimax,BGMjmax,BGMkmax)) &
                     ,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/FIBGM_offsetElem_Shared(BGMimax,BGMjmax,BGMkmax)+FIBGM_nElem_Shared(BGMimax,BGMjmax,BGMkmax)/)&
                     ,FIBGM_Element_Shared_Win,FIBGM_Element_Shared)
CALL MPI_WIN_LOCK_ALL(0,FIBGM_Element_Shared_Win,IERROR)

DO kBGM = BGMkmin,BGMkmax
  DO jBGM = BGMjmin,BGMjmax
    DO iBGM = BGMimin,BGMimax
      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

DO iElem = firstHaloElem, lastHaloElem
  ElemID = offsetHaloElem(iElem)
  IF (ElemInfo_Shared(ElemInfoSize+1,ElemID).EQ.0) CYCLE
  BGMCellXmin = ElemToBGM(1,ElemID)
  BGMCellXmax = ElemToBGM(2,ElemID)
  BGMCellYmin = ElemToBGM(3,ElemID)
  BGMCellYmax = ElemToBGM(4,ElemID)
  BGMCellZmin = ElemToBGM(5,ElemID)
  BGMCellZmax = ElemToBGM(6,ElemID)
  ! add current Element to BGM-Elem
  DO kBGM = BGMCellZmin,BGMCellZmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO iBGM = BGMCellXmin,BGMCellXmax
        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
        FIBGM_Element_Shared( FIBGM_offsetElem_Shared(iBGM,jBGM,kBGM) & ! offset of BGM cell in 1D array
                            + offsetElemsInBGMCell(iBGM,jBGM,kBGM)    & ! offset of BGM nElems in local proc
                            + GEO%FIBGM(iBGM,jBGM,kBGM)%nElem         ) = ElemID
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem = firstHaloElem, lastHaloElem
DO iElem = offsetElem+1, offsetElem+nElems
  BGMCellXmin = ElemToBGM(1,iElem)
  BGMCellXmax = ElemToBGM(2,iElem)
  BGMCellYmin = ElemToBGM(3,iElem)
  BGMCellYmax = ElemToBGM(4,iElem)
  BGMCellZmin = ElemToBGM(5,iElem)
  BGMCellZmax = ElemToBGM(6,iElem)
  ! add current Element to BGM-Elem
  DO kBGM = BGMCellZmin,BGMCellZmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO iBGM = BGMCellXmin,BGMCellXmax
        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
        FIBGM_Element_Shared( FIBGM_offsetElem_Shared(iBGM,jBGM,kBGM) & ! offset of BGM cell in 1D array
                            + offsetElemsInBGMCell(iBGM,jBGM,kBGM)    & ! offset of BGM nElems in local proc
                            + GEO%FIBGM(iBGM,jBGM,kBGM)%nElem         ) = iElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem

!--- map elements to background cells
! alternative if nElem is counted with cycles
!DO iElem = firstElem, lastElem
!  IF (ElemInfo_Shared(ElemInfoSize+1,iElem).EQ.0) CYCLE
!  BGMCellXmin = ElemToBGM(1,iElem)
!  BGMCellXmax = ElemToBGM(2,iElem)
!  BGMCellYmin = ElemToBGM(3,iElem)
!  BGMCellYmax = ElemToBGM(4,iElem)
!  BGMCellZmin = ElemToBGM(5,iElem)
!  BGMCellZmax = ElemToBGM(6,iElem)
!  ! add current Element to BGM-Elem
!  DO kBGM = BGMCellZmin,BGMCellZmax
!    DO jBGM = BGMCellYmin,BGMCellYmax
!      DO iBGM = BGMCellXmin,BGMCellXmax
!        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
!        FIBGM_Element_Shared( FIBGM_offsetElem_Shared(iBGM,jBGM,kBGM) & ! offset of BGM cell in 1D array
!                            + offsetElemsInBGMCell(iBGM,jBGM,kBGM)    & ! offset of BGM nElems in local proc
!                            + GEO%FIBGM(iBGM,jBGM,kBGM)%nElem         ) = iElem
!      END DO ! kBGM
!    END DO ! jBGM
!  END DO ! iBGM
!END DO ! iElem
DEALLOCATE(offsetElemsInBGMCell)

CALL MPI_WIN_SYNC(FIBGM_Element_Shared,IERROR)
CALL MPI_BARRIER(MPI_COMM_SHARED,iError)

#else  /*NOT USE_MPI*/

DO iElem = 1, nElems
  xmin=MINVAL(NodeCoords(1,:,:,:,iElem))
  xmax=MAXVAL(NodeCoords(1,:,:,:,iElem))
  ymin=MINVAL(NodeCoords(2,:,:,:,iElem))
  ymax=MAXVAL(NodeCoords(2,:,:,:,iElem))
  zmin=MINVAL(NodeCoords(3,:,:,:,iElem))
  zmax=MAXVAL(NodeCoords(3,:,:,:,iElem))
  ElemToBGM(1,iElem) = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
  ElemToBGM(2,iElem) = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
  ElemToBGM(3,iElem) = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
  ElemToBGM(4,iElem) = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
  ElemToBGM(5,iElem) = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
  ElemToBGM(6,iElem) = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
END DO ! iElem = 1, nElems

DO kBGM = BGMkmin,BGMkmax
  DO jBGM = BGMjmin,BGMjmax
    DO iBGM = BGMimin,BGMimax
      IF(GEO%FIBGM(iBGM,jBGM,kBGM)%nElem.EQ.0) CYCLE
      ALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element(1:GEO%FIBGM(iBGM,jBGM,kBGM)%nElem))
      GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
    END DO ! kBGM
  END DO ! jBGM
END DO ! iBGM

!--- map elements to background cells
DO iElem=1,PP_nElems
  BGMCellXmin = ElemToBGM(1,iElem)
  BGMCellXmax = ElemToBGM(2,iElem)
  BGMCellYmin = ElemToBGM(3,iElem)
  BGMCellYmax = ElemToBGM(4,iElem)
  BGMCellZmin = ElemToBGM(5,iElem)
  BGMCellZmax = ElemToBGM(6,iElem)
  ! add current Element to BGM-Elem
  DO kBGM = BGMCellZmin,BGMCellZmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO iBGM = BGMCellXmin,BGMCellXmax
        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
        GEO%FIBGM(iBGM,jBGM,kBGM)%Element(GEO%FIBGM(iBGM,jBGM,kBGM)%nElem) = iElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem
DEALLOCATE(NodeCoords)
#endif  /*USE_MPI*/

END SUBROUTINE BuildBGM


END MODULE MOD_Particle_BGM
