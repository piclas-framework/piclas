#include "boltzplatz.h"

MODULE MOD_Liquid_Boundary
!===================================================================================================================================
! defines routines and functions for liquid boundaries
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
INTERFACE Init_Liquid_Boundary
  MODULE PROCEDURE Init_Liquid_Boundary
END INTERFACE
INTERFACE Evaporation
  MODULE PROCEDURE Evaporation
END INTERFACE

PUBLIC :: Evaporation,Init_Liquid_Boundary
!===================================================================================================================================

CONTAINS

SUBROUTINE Init_Liquid_Boundary()
!===================================================================================================================================
! init of particle liquid boundary interfaces
! 1) mark sides for liquid boundary consideration
! 2) init MPI communication
!===================================================================================================================================
USE MOD_Particle_Vars,          ONLY : nSpecies
USE MOD_DSMC_Vars,              ONLY : Liquid
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   INTEGER                          :: iSpec
   REAL                             :: PartAds, PartEvap, RanNum, Tpois
!===================================================================================================================================
! allocate info and constants
ALLOCATE( Liquid%Info(1:nSpecies))
! initialize info and constants
DO iSpec = 1,nSpecies
#if (PP_TimeDiscMethod==42)
  Liquid%Info(iSpec)%MeanProbAds = 0.
  Liquid%Info(iSpec)%MeanProbDes = 0.
  Liquid%Info(iSpec)%MeanEads = 0.
  Liquid%Info(iSpec)%WallCollCount = 0
  Liquid%Info(iSpec)%NumOfAds = 0
  Liquid%Info(iSpec)%NumOfDes = 0
#endif
END DO
! allocate and initialize liquid variables
ALLOCATE( Liquid%ProbCondens(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Liquid%ProbEvap(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Liquid%SumCondensPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nTotalSides,1:nSpecies),&
          Liquid%SumEvapPart(1:nSurfSample,1:nSurfSample,1:SurfMesh%nSides,1:nSpecies))
Liquid%ProbCondens(:,:,:,:) = 1. !0.
Liquid%ProbEvap(:,:,:,:) = 0.
Liquid%SumCondensPart(:,:,:,:) = 0
Liquid%SumEvapPart(:,:,:,:) = 0

END SUBROUTINE Init_Liquid_Boundary


SUBROUTINE Evaporation()
!===================================================================================================================================
! calculation of evaporating particle number when particles are deleted at condensation and inserted at evaporation
!===================================================================================================================================
USE MOD_Particle_Vars,          ONLY : nSpecies, BoltzmannConst, Species
USE MOD_DSMC_Vars,              ONLY : Liquid
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh, PartBound
USE MOD_Mesh_Vars,              ONLY : BC
USE MOD_TimeDisc_Vars,          ONLY : dt
!===================================================================================================================================
  IMPLICIT NONE
!===================================================================================================================================
! Local variable declaration
   REAL, PARAMETER                  :: PI=3.14159265358979323846_8
   INTEGER                          :: iSurfSide, iSpec, p, q, Npois
   REAL                             :: PartEvap, RanNum, Tpois
   REAL                             :: LiquidTemp
   REAL                             :: pressure_vapor, A, B, C
!===================================================================================================================================
#if (PP_TimeDiscMethod==42)
  Liquid%Info(:)%NumOfDes = 0
#endif
  Liquid%SumEvapPart(:,:,:,:) = 0
  DO iSpec = 1,nSpecies
    DO iSurfSide = 1,SurfMesh%nSides
      IF (PartBound%SolidState(PartBound%MapToPartBC(BC( SurfMesh%SurfSideToGlobSideMap(iSurfSide) )))) CYCLE
      IF (PartBound%LiquidSpec(PartBound%MapToPartBC(BC( SurfMesh%SurfSideToGlobSideMap(iSurfSide) ))).NE.iSpec) CYCLE
      LiquidTemp = PartBound%WallTemp(PartBound%MapToPartBC(BC(SurfMesh%SurfSideToGlobSideMap(iSurfSide))))
      DO q = 1,nSurfSample
        DO p = 1,nSurfSample
          A = 5.40221
          B = 1838.675
          C = -31.737
          pressure_vapor = 10 ** (A- B/(C+LiquidTemp)) * 1e5 !transformation bar -> Pa
          PartEvap = pressure_vapor / ( 2*PI*Species(iSpec)%MassIC*BoltzmannConst*LiquidTemp)**0.5 &
                   * SurfMesh%SurfaceArea(p,q,iSurfSide) / Species(iSpec)%MacroParticleFactor * dt
          CALL RANDOM_NUMBER(RanNum)            
          IF (EXP(-PartEvap).LE.TINY(PartEvap)) THEN
            Liquid%SumEvapPart(p,q,iSurfSide,iSpec) = Liquid%SumEvapPart(p,q,iSurfSide,iSpec)&
            +INT(PartEvap + RanNum)
          ELSE !poisson-sampling instead of random rounding (reduces numeric non-equlibrium effects [Tysanner and Garcia 2004]
            Npois=0
            Tpois=1.0
            DO
              Tpois=RanNum*Tpois
              IF (Tpois.LT.TINY(Tpois)) THEN
                Liquid%SumEvapPart(p,q,iSurfSide,iSpec) = Liquid%SumEvapPart(p,q,iSurfSide,iSpec)&
                +INT(PartEvap + RanNum)
                EXIT
              END IF
              IF (Tpois.GT.EXP(-PartEvap)) THEN
                Npois=Npois+1
                CALL RANDOM_NUMBER(RanNum)
              ELSE
                Liquid%SumEvapPart(p,q,iSurfSide,iSpec) = Liquid%SumEvapPart(p,q,iSurfSide,iSpec)+Npois
                EXIT
              END IF
            END DO
          END IF
#if (PP_TimeDiscMethod==42)
          Liquid%Info(iSpec)%NumOfDes = Liquid%Info(iSpec)%NumOfDes + Liquid%SumEvapPart(p,q,iSurfSide,iSpec)
#endif
        END DO  
      END DO
    END DO
  END DO
  
END SUBROUTINE Evaporation

#ifdef MPI
SUBROUTINE ExchangeCondensNum() 
!===================================================================================================================================
! exchange the number of condensing particles on halo surface 
! only processes with samling sides in their halo region and the original process participate on the communication
! structure is similar to surface sampling/particle communication
! each process sends his halo-information directly to the origin process by use of a list, containing the surfsideids for sending
! the receiving process adds the new data to his own sides
!===================================================================================================================================
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_Particle_Vars               ,ONLY:nSpecies
USE MOD_DSMC_Vars                   ,ONLY:Liquid
USE MOD_Particle_Boundary_Vars      ,ONLY:SurfComm,nSurfSample
USE MOD_Particle_MPI_Vars           ,ONLY:CondensSendBuf,CondensRecvBuf,SurfExchange
!----------------------------------------------------------------------------------------------------------------------------------!
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
! INPUT VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------!
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: MessageSize,nValues,iSurfSide,SurfSideID
INTEGER                         :: iPos,p,q,iProc
INTEGER                         :: recv_status_list(1:MPI_STATUS_SIZE,1:SurfCOMM%nMPINeighbors)
!===================================================================================================================================

nValues=nSpecies*(nSurfSample)**2

! open receive buffer
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesRecv(iProc)*nValues
  CALL MPI_IRECV( CondensRecvBuf(iProc)%content_int             &
                , MessageSize                                  &
                , MPI_INT                                      &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID     &
                , 1011                                         &
                , SurfCOMM%COMM                                &
                , SurfExchange%RecvRequest(iProc)              & 
                , IERROR )
END DO ! iProc

! build message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  iPos=0
  DO iSurfSide=1,SurfExchange%nSidesSend(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%SendList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        CondensSendBuf(iProc)%content_int(iPos+1:iPos+nSpecies)= Liquid%SumCondensPart(p,q,SurfSideID,:)
        iPos=iPos+nSpecies
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
!     Adsorption%SumAdsorbPart(:,:,SurfSideID,:)=0.
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
END DO

! send message
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
  CALL MPI_ISEND( CondensSendBuf(iProc)%content_int         &
                , MessageSize                              & 
                , MPI_INT                                  &
                , SurfCOMM%MPINeighbor(iProc)%NativeProcID & 
                , 1011                                     &
                , SurfCOMM%COMM                            &   
                , SurfExchange%SendRequest(iProc)          &
                , IERROR )                                     
END DO ! iProc                                                

! 4) Finish Received number of particles
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesSend(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%SendRequest(iProc),MPIStatus,IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
  IF(SurfExchange%nSidesRecv(iProc).NE.0) THEN
    CALL MPI_WAIT(SurfExchange%RecvRequest(iProc),recv_status_list(:,iProc),IERROR)
    IF(IERROR.NE.MPI_SUCCESS) CALL abort(&
__STAMP__&
          ,' MPI Communication error', IERROR)
  END IF
END DO ! iProc

! add data do my list
DO iProc=1,SurfCOMM%nMPINeighbors
  IF(SurfExchange%nSidesRecv(iProc).EQ.0) CYCLE
  MessageSize=SurfExchange%nSidesSend(iProc)*nValues
  iPos=0
  DO iSurfSide=1,SurfExchange%nSidesRecv(iProc)
    SurfSideID=SurfCOMM%MPINeighbor(iProc)%RecvList(iSurfSide)
    DO q=1,nSurfSample
      DO p=1,nSurfSample
        Liquid%SumCondensPart(p,q,SurfSideID,:)=Liquid%SumCondensPart(p,q,SurfSideID,:) &
                                         +CondensRecvBuf(iProc)%content_int(iPos+1:iPos+nSpecies)
        iPos=iPos+nSpecies
      END DO ! p=0,nSurfSample
    END DO ! q=0,nSurfSample
  END DO ! iSurfSide=1,nSurfExchange%nSidesSend(iProc)
  CondensRecvBuf(iProc)%content_int = 0
  CondensSendBuf(iProc)%content_int = 0
END DO ! iProc

END SUBROUTINE ExchangeCondensNum
#endif /*MPI*/

SUBROUTINE Finalize_Liquid_Boundary()
!===================================================================================================================================
! Deallocate liquid boundary vars
!===================================================================================================================================
! MODULES
USE MOD_DSMC_Vars,              ONLY : Liquid
USE MOD_Particle_Boundary_Vars, ONLY : nSurfSample, SurfMesh
! #ifdef MPI
! USE MOD_Particle_Boundary_Vars, ONLY : SurfCOMM
! USE MOD_Particle_MPI_Vars,      ONLY : SurfExchange
! USE MOD_Particle_MPI_Vars,      ONLY : AdsorbSendBuf,AdsorbRecvBuf,SurfDistSendBuf,SurfDistRecvBuf
! #endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
! INTEGER                      :: subsurfxi,subsurfeta,iSurfSide
! #ifdef MPI
! INTEGER                      :: iProc
! #endif /*MPI*/
!===================================================================================================================================
SDEALLOCATE(Liquid%Info)

SDEALLOCATE(Liquid%ProbCondens)
SDEALLOCATE(Liquid%ProbEvap)
SDEALLOCATE(Liquid%SumCondensPart)
SDEALLOCATE(Liquid%SumEvapPart)

END SUBROUTINE Finalize_Liquid_Boundary

END MODULE MOD_Liquid_Boundary
