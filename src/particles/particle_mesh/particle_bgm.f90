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
!-----------------------------------------------------------------------------------------------------------------------------------
! required variables
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

PUBLIC::CountPartsPerElem
!===================================================================================================================================
!
PUBLIC::DefineParametersParticleMesh
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
! computes the element indices of an given element in the BGM-mesh
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
#if USE_MPI
USE MOD_Mesh_Vars            ,ONLY: offsetElem,nElems
USE MOD_MPI_Vars             ,ONLY: nMPISides_Proc,nNbProcs,NbProc
USE MOD_MPI_Shared_Vars
USE MOD_MPI_Shared           ,ONLY: Allocate_Shared
#else
USE MOD_Mesh_Vars            ,ONLY: NodeCoords, nElems
#endif /*USE_MPI*/
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: ElemID
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem
INTEGER,ALLOCATABLE :: ElemToBGM(:,:)
!===================================================================================================================================

ALLOCATE(ElemToBGM(1:6,nElems))
DO iElem = 1, nElems
#if USE_MPI
  firstNodeID=ElemInfo_Shared(ELEM_FirstNodeInd,offsetElem+iElem)
  nNodeIDs=ElemInfo_Shared(ELEM_LastNodeInd,offsetElem+iElem)-ElemInfo_Shared(ELEM_FirstNodeInd,offsetElem+iElem)

  xmin=MINVAL(NodeCoords_Shared(1,firstNodeID:firstNodeID+nNodeIDs))
  xmax=MAXVAL(NodeCoords_Shared(1,firstNodeID:firstNodeID+nNodeIDs))
  ymin=MINVAL(NodeCoords_Shared(2,firstNodeID:firstNodeID+nNodeIDs))
  ymax=MAXVAL(NodeCoords_Shared(2,firstNodeID:firstNodeID+nNodeIDs))
  zmin=MINVAL(NodeCoords_Shared(3,firstNodeID:firstNodeID+nNodeIDs))
  zmax=MAXVAL(NodeCoords_Shared(3,firstNodeID:firstNodeID+nNodeIDs))
#else
  xmin=MINVAL(NodeCoords(1,:,:,:,iElem))
  xmax=MAXVAL(NodeCoords(1,:,:,:,iElem))
  ymin=MINVAL(NodeCoords(2,:,:,:,iElem))
  ymax=MAXVAL(NodeCoords(2,:,:,:,iElem))
  zmin=MINVAL(NodeCoords(3,:,:,:,iElem))
  zmax=MAXVAL(NodeCoords(3,:,:,:,iElem))
#endif /*USE_MPI*/
  ElemToBGM(1,iElem) = CEILING((xmin-GEO%xminglob)/GEO%FIBGMdeltas(1))
  ElemToBGM(2,iElem) = CEILING((xmax-GEO%xminglob)/GEO%FIBGMdeltas(1))
  ElemToBGM(3,iElem) = CEILING((ymin-GEO%yminglob)/GEO%FIBGMdeltas(2))
  ElemToBGM(4,iElem) = CEILING((ymax-GEO%yminglob)/GEO%FIBGMdeltas(2))
  ElemToBGM(5,iElem) = CEILING((zmin-GEO%zminglob)/GEO%FIBGMdeltas(3))
  ElemToBGM(6,iElem) = CEILING((zmax-GEO%zminglob)/GEO%FIBGMdeltas(3))
END DO ! iElem = 1, nElems

#if USE_MPI
MPISharedSize = INT(6*nTotalElems,MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/6,nTotalElems/),ElemToBGM_Shared_Win,ElemToBGM_Shared)
CALL MPI_WIN_LOCK_ALL(0,ElemToBGM_Shared_Win,IERROR)
ElemToBGM_Shared(:,offsetElem:offsetElem+nElems) = ElemToBGM(:,:)
CALL MPI_WIN_SYNC(ElemToBGM_Shared_Win,IERROR)
IF(myRank_Shared.EQ.0)THEN
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ElemToBGM_Shared,6*nTotalElems,MPI_INTEGER,MPI_COMM_LEADERS_SHARED,IERROR)
END IF
#endif /*USE_MPI*/



CALL InitPeriodicBC()

! deallocate stuff // required for dynamic load balance
#if USE_MPI
IF (ALLOCATED(GEO%FIBGM)) THEN
  DO iBGM=GEO%FIBGMimin,GEO%FIBGMimax
    DO jBGM=GEO%FIBGMjmin,GEO%FIBGMjmax
      DO kBGM=GEO%FIBGMkmin,GEO%FIBGMkmax
        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%Element)
        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%ShapeProcs)
        SDEALLOCATE(GEO%FIBGM(iBGM,jBGM,kBGM)%PaddingProcs)
!           SDEALLOCATE(GEO%FIBGM(i,k,l)%SharedProcs)
      END DO
    END DO
  END DO
  DEALLOCATE(GEO%FIBGM)
END IF
#endif /*USE_MPI*/

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
    halo_eps = haloe_eps + r_sf
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

#if USE_MPI
! enlarge the BGM grid for safety reasons
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

GEO%FIBGMimax=BGMimax
GEO%FIBGMimin=BGMimin
GEO%FIBGMjmax=BGMjmax
GEO%FIBGMjmin=BGMjmin
GEO%FIBGMkmax=BGMkmax
GEO%FIBGMkmin=BGMkmin

ALLOCATE(GEO%FIBGM(BGMimin:BGMimax,BGMjmin:BGMjmax,BGMkmin:BGMkmax))

! null number of element per BGM cell
DO kBGM = BGMkmin,BGMkmax
   DO jBGM = BGMjmin,BGMjmax
     DO iBGM = BGMimin,BGMimax
         GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = 0
      END DO
   END DO
END DO

!--- compute number of elements in each background cell
DO iElem=1,PP_nElems
  ! here fancy stuff, because element could be wide out of element range
  BGMCellXmin = ElemToBGM_Shared(1,offsetElem+iElem)
  BGMCellXmax = ElemToBGM_Shared(2,offsetElem+iElem)
  BGMCellYmin = ElemToBGM_Shared(3,offsetElem+iElem)
  BGMCellYmax = ElemToBGM_Shared(4,offsetElem+iElem)
  BGMCellZmin = ElemToBGM_Shared(5,offsetElem+iElem)
  BGMCellZmax = ElemToBGM_Shared(6,offsetElem+iElem)
  ! add current element to number of BGM-elems
  DO iBGM = BGMCellXmin,BGMCellXmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO kBGM = BGMCellZmin,BGMCellZmax
         GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem

! find max nelems
CALL MPI_EXSCAN(...)
! bcast of total bgm elems
CALL MPI_BCAST(bgm%totalelems
#if USE_MPI
MPISharedSize = INT((BGMimax-BGMimin)*(BGMjmax-BGMjmin)*(BGMkmax-BGMkmin),MPI_ADDRESS_KIND)*MPI_ADDRESS_KIND
CALL Allocate_Shared(MPISharedSize,(/BGMimax-BGMimin,BGMjmax-BGMjmin,BGMkmax-BGMkmin/),GEO%FIBGM_Shared_Win,GEO%FIBGM_Shared)
CALL MPI_WIN_LOCK_ALL(0,GEO%FIBGM_Shared_Win,IERROR)
#endif  /*USE_MPI*/


!--- allocate mapping variable and clean number for mapping (below)
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
  ! here fancy stuff, because element could be wide out of element range
  BGMCellXmin = ElemToBGM_Shared(1,offsetElem+iElem)
  BGMCellXmax = ElemToBGM_Shared(2,offsetElem+iElem)
  BGMCellYmin = ElemToBGM_Shared(3,offsetElem+iElem)
  BGMCellYmax = ElemToBGM_Shared(4,offsetElem+iElem)
  BGMCellZmin = ElemToBGM_Shared(5,offsetElem+iElem)
  BGMCellZmax = ElemToBGM_Shared(6,offsetElem+iElem)
  ! add current Element to BGM-Elem
  DO kBGM = BGMCellZmin,BGMCellZmax
    DO jBGM = BGMCellYmin,BGMCellYmax
      DO iBGM = BGMCellXmin,BGMCellXmax
        GEO%FIBGM(iBGM,jBGM,kBGM)%nElem = GEO%FIBGM(iBGM,jBGM,kBGM)%nElem + 1
        GEO%FIBGM(iBGM,jBGM,kBGM)%Element(offsetElem+GEO%FIBGM(iBGM,jBGM,kBGM)%nElem) = iElem
      END DO ! kBGM
    END DO ! jBGM
  END DO ! iBGM
END DO ! iElem

END SUBROUTINE BuildBGM
