#include "boltzplatz.h"

MODULE MOD_Particle_Boundary_Sampling
!===================================================================================================================================
!! Determines how particles interact with a given boundary condition 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitParticleBoundarySampling
  MODULE PROCEDURE InitParticleBoundarySampling
END INTERFACE

INTERFACE FinalizeParticleBoundarySampling
  MODULE PROCEDURE FinalizeParticleBoundarySampling
END INTERFACE


PUBLIC::InitParticleBoundarySampling
PUBLIC::FinalizeParticleBoundarySampling
!===================================================================================================================================

CONTAINS

SUBROUTINE InitParticleBoundarySampling() 
!===================================================================================================================================
! init of particle boundary sampling
! default: use for sampling same polynomial degree as NGeo
! 1) mark sides for sampling
! 2) build special MPI communicator
#ifdef MPI
#endif /*MPI*/
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Mesh_Vars               ,ONLY:NGeo,BC
USE MOD_ReadInTools             ,ONLY:GETINT
USE MOD_Particle_Boundary_Vars  ,ONLY:nSurfSample,dXiEQ_SurfSample,PartBound,XiEQ_SurfSample,SurfMesh
USE MOD_Particle_Mesh_Vars      ,ONLY:nTotalSides
USE MOD_MPI_Vars                ,ONLY:offsetSurfElemMPI
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: q
CHARACTER(2)                :: hilf
!===================================================================================================================================
 
SWRITE(UNIT_stdOut,'(A)') ' INIT SURFACE SAMPLING ...'
WRITE(UNIT=hilf,FMT='(I2)') NGeo
nSurfSample = GETINT('DSMC-nSurfSample',hilf)
 
ALLOCATE(XiEQ_SurfSample(0:nSurfSample))

dXiEQ_SurfSample =2./REAL(nSurfSample)
DO q=0,nSurfSample
  XiEQ_SurfSample(q) = dXiEQ_SurfSample * REAL(q) - 1. 
END DO

! change basis or recompute???

! first, get number of bc sides
!SurfMesh%nSurfaceBCSides = 0
!DO iSide=1,nTotalSides
!  IF(BC(iSide).LE.1) CYCLE
!  IF (PartBound%TargetBoundCond(PartBound%MapToPartBC(BC(iSide))).EQ.PartBound%ReflectiveBC) THEN  
!    SurfMesh%nSurfaceBCSides = SurfMesh%nSurfaceBCSides + 1
!    SurfMesh%SideToSurfID(iSide) = SurfMesh%nSurfaceBCSides
!  END IF
!END DO ! iSide=1,nTotalSides



! calculating offset of surface elements for DSMC surface output

!#ifdef MPI
!
!SDEALLOCATE(offsetSurfElemMPI)
!ALLOCATE(offsetSurfElemMPI(0:nProcessors))
!offsetSurfElemMPI=0
!
!countSurfElem=0
!
!DO iSide=1,nBCSides
!  IF (BoundaryType(BC(iSide),1).EQ.4) THEN
!    countSurfElem = countSurfElem + 1
!  END IF
!END DO
!
!!IF (MPIroot) THEN
!ALLOCATE(countSurfElemMPI(0:nProcessors-1))
!countSurfElemMPI=0
!!END IF
!
!CALL MPI_GATHER(countSurfElem,1,MPI_INTEGER,countSurfElemMPI,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
!
!IF (MPIroot) THEN
!DO iProc=1,nProcessors-1
!  offsetSurfElemMPI(iProc)=SUM(countSurfElemMPI(0:iProc-1))
!END DO
!offsetSurfElemMPI(nProcessors)=SUM(countSurfElemMPI(:))
!END IF
!
!CALL MPI_BCAST (offsetSurfElemMPI,size(offsetSurfElemMPI),MPI_INTEGER,0,MPI_COMM_WORLD,iError)
!
!offsetSurfElem=offsetSurfElemMPI(myRank)
!
!DEALLOCATE(countSurfElemMPI)
!#else /* MPI */
!offsetSurfElem=0          ! offset is the index of first entry, hdf5 array starts at 0-.GT. -1 
!#endif /* MPI */





#ifdef MPI



#endif /*MPI*/

SWRITE(UNIT_stdOut,'(A)') ' ... DONE.'

END SUBROUTINE InitParticleBoundarySampling


SUBROUTINE FinalizeParticleBoundarySampling() 
!===================================================================================================================================
! deallocate everything
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Particle_Boundary_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(XiEQ_SurfSample)

END SUBROUTINE FinalizeParticleBoundarySampling


END MODULE MOD_Particle_Boundary_Sampling
