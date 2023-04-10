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

MODULE MOD_Restart_Tools
!===================================================================================================================================
! Module containing tools/procedures for handling PICLas restarts
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------

#if USE_HDG
PUBLIC :: RecomputeLambda
#endif /*USE_HDG*/
!===================================================================================================================================

CONTAINS


#if USE_HDG
!===================================================================================================================================
!> The lambda-solution is stored per side, however, the side-list is computed with the OLD domain-decomposition. To allow for
!> a change in the load-distribution, number of used cores, etc,... lambda has to be recomputed ONCE
!===================================================================================================================================
SUBROUTINE RecomputeLambda(t)
! MODULES
use MOD_Globals
USE MOD_DG_Vars            ,ONLY: U
USE MOD_PreProc
USE MOD_HDG                ,ONLY: HDG
USE MOD_TimeDisc_Vars      ,ONLY: iter
#ifdef PARTICLES
USE MOD_PICDepo            ,ONLY: Deposition
USE MOD_HDG_Vars           ,ONLY: UseBRElectronFluid,BRElectronsRemoved
#if USE_MPI
USE MOD_Particle_MPI       ,ONLY: IRecvNbOfParticles, MPIParticleSend,MPIParticleRecv,SendNbOfparticles
#endif /*USE_MPI*/
#endif /*PARTICLES*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

#ifdef PARTICLES
! Deposition of particles
CALL Deposition()
#endif /*PARTICLES*/

! recompute fields
! EM field
#ifdef PARTICLES
IF(UseBRElectronFluid.AND.BRElectronsRemoved)THEN
  ! When using BR electron fluid model, all electrons are removed from the restart file
  CALL HDG(t,U,iter,ForceCGSolverIteration_opt=.TRUE.)
ELSE
#endif /*PARTICLES*/
  CALL HDG(t,U,iter)
#ifdef PARTICLES
END IF ! UseBRElectronFluid
#endif /*PARTICLES*/

END SUBROUTINE RecomputeLambda
#endif /*USE_HDG*/


END MODULE MOD_Restart_Tools
