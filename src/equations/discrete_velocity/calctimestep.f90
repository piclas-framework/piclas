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

MODULE MOD_CalcTimeStep
!===================================================================================================================================
! Low-Storage Runge-Kutta integration of degree 3 for one step.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE CALCTIMESTEP
  MODULE PROCEDURE CALCTIMESTEP
END INTERFACE


PUBLIC :: CALCTIMESTEP
!===================================================================================================================================

CONTAINS

FUNCTION CALCTIMESTEP()
!===================================================================================================================================
! Calculate the time step for the current update of U for DVM
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Particle_Mesh_Vars,ONLY: ElemVolume_Shared
USE MOD_Mesh_Vars_FV      ,ONLY: NormVec_FV, SurfElem_FV
USE MOD_Mesh_Vars         ,ONLY: ElemToSide
USE MOD_Equation_Vars_FV  ,ONLY: DVMSpecData, DVMnSpecies
USE MOD_TimeDisc_Vars     ,ONLY: CFLScale
#ifdef PARTICLES
USE MOD_Mesh_Vars         ,ONLY: offsetElem
USE MOD_Mesh_Tools        ,ONLY: GetCNElemID
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                         :: CalcTimeStep
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iElem,CNELemID,locSideID,SideID,flip
INTEGER                      :: iVel,jVel,kVel,upos,iSpec,vFirstID
REAL                         :: FluxFac(PP_nVar_FV), SideFac
REAL                         :: TimeStepConv,locTimeStepConv
REAL                         :: n_loc(3)
!===================================================================================================================================
locTimeStepConv=HUGE(1.)
DO iElem=1,PP_nElems
#if USE_MPI && defined(PARTICLES)
  CNElemID=GetCNElemID(iElem+offSetElem)
#else
  CNElemID=iElem
#endif
  FluxFac = 0.
#if PP_dim == 3
  DO locSideID=1,6
#else
  DO locSideID=2,5
#endif
    SideID=ElemToSide(E2S_SIDE_ID,locSideID,iElem)
    flip=ElemToSide(E2S_FLIP,locSideID,iElem)
    IF (flip.EQ.0) THEN !master side, outgoing normal OK
      n_loc = NormVec_FV(:,0,0,SideID)
    ELSE !slave side, flip ingoing normal
      n_loc = -NormVec_FV(:,0,0,SideID)
    END IF
    SideFac = SurfElem_FV(0,0,SideID)/ElemVolume_Shared(CNElemID)
    vFirstID = 0
    DO iSpec =1,DVMnSpecies
      DO kVel=1, DVMSpecData(iSpec)%nVelos(3);   DO jVel=1, DVMSpecData(iSpec)%nVelos(2);   DO iVel=1, DVMSpecData(iSpec)%nVelos(1)
        upos=iVel+(jVel-1)*DVMSpecData(iSpec)%nVelos(1)+(kVel-1)*DVMSpecData(iSpec)%nVelos(1)*DVMSpecData(iSpec)%nVelos(2)+vFirstID
        FluxFac(upos) = FluxFac(upos) + SideFac*MAX(0.,n_loc(1)*DVMSpecData(iSpec)%Velos(iVel,1) &
                                                      +n_loc(2)*DVMSpecData(iSpec)%Velos(jVel,2) &
                                                      +n_loc(3)*DVMSpecData(iSpec)%Velos(kVel,3))
      END DO; END DO; END DO
      vFirstID = vFirstID + DVMSpecData(iSpec)%nVar
    END DO
  END DO
  locTimeStepConv = MIN(locTimeStepConv,CFLScale/MAXVAL(FluxFac))

  IF(locTimeStepConv.NE.locTimeStepConv)THEN
    ERRWRITE(*,'(A,I0,A,I0,A,I0,A,I0)')'Convective timestep NaN on proc ',myRank,' at global position (iElem): ',iElem
    ERRWRITE(*,*)'dt_conv=',locTimeStepConv
    CALL abort(&
        __STAMP__&
        ,'Convective timestep NaN!',999,999.)
  END IF
END DO ! iElem
#if USE_MPI
CALL MPI_ALLREDUCE(locTimeStepConv,TimeStepConv,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_PICLAS,iError)
#else
TimeStepConv=locTimeStepConv
#endif /*USE_MPI*/
CalcTimeStep=TimeStepConv

END FUNCTION CALCTIMESTEP

END MODULE MOD_CalcTimeStep
