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
#if USE_QDS_DG
#include "piclas.h"

MODULE MOD_QDS_DG
!===================================================================================================================================
!> Contains the routines to
!> - initialize the QDS-DG method
!> - finalize the QDS-DG method variables
!> - determine the weak form time derivative for the QDS-G method
!> - calculate the QDS-MacroValues from the QDS-DG variables (QDSReCalculateDGValues)
!> - calculate the QDS-DG variables form the QDS-MacroValues (QDSCalculateMacroValues)
!===================================================================================================================================
! MODULES
!USE MOD_io_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE QDS_InitDG
  MODULE PROCEDURE QDS_InitDG
END INTERFACE
INTERFACE QDSTimeDerivative
  MODULE PROCEDURE QDSTimeDerivative
END INTERFACE
INTERFACE QDS_FinalizeDG
  MODULE PROCEDURE QDS_FinalizeDG
END INTERFACE
INTERFACE QDSReCalculateDGValues
  MODULE PROCEDURE QDSReCalculateDGValues
END INTERFACE
INTERFACE QDSCalculateMacroValues
  MODULE PROCEDURE QDSCalculateMacroValues
END INTERFACE

PUBLIC::QDS_InitDG
PUBLIC::QDS_FinalizeDG
PUBLIC::QDSTimeDerivative
PUBLIC::QDSReCalculateDGValues
PUBLIC::QDSCalculateMacroValues
!===================================================================================================================================
CONTAINS


SUBROUTINE QDS_InitDG
!===================================================================================================================================
!> Allocate all QDS variables, determine
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_QDS_DG_Vars
USE MOD_Globals_Vars,       ONLY:BoltzmannConst
USE MOD_Globals,            ONLY:abort,UNIT_stdOut,mpiroot
USE MOD_ReadInTools,        ONLY:GETLOGICAL
USE MOD_Mesh_Vars,          ONLY:nSides
USE MOD_Globals_Vars,       ONLY:PI
USE MOD_Restart_Vars,       ONLY:DoRestart
USE MOD_QDS_Equation_vars,  ONLY:QDSnVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: tempNorm
INTEGER :: iWeight
REAL    :: Velo(3), Temp, Dens, Mass
!===================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT QDS-DG...'
IF(QDSInitDGIsDone)THEN
  CALL abort(&
      __STAMP__&
      ,'QDS_InitDG (DG method) not ready to be called or already called.')
END IF
DoQDS                      = GETLOGICAL('DoQDS','.FALSE.')
IF(DoQDS)THEN
  QDSnVar_macro=6
  nQDSElems=PP_nElems

  ALLOCATE(UQDS       (1:QDSnVar,0:PP_N,0:PP_N,0:PP_N,1:nQDSElems))
  ALLOCATE(UQDSt      (1:QDSnVar,0:PP_N,0:PP_N,0:PP_N,1:nQDSElems))
  UQDS=0.
  UQDSt=0.

  ALLOCATE(UQDS_master(1:QDSnVar,0:PP_N,0:PP_N,1:nSides))
  ALLOCATE(UQDS_slave( 1:QDSnVar,0:PP_N,0:PP_N,1:nSides))
  UQDS_master=0.
  UQDS_slave=0.

  ALLOCATE(FluxQDS_Master(1:QDSnVar,0:PP_N,0:PP_N,1:nSides))
  ALLOCATE(FluxQDS_Slave( 1:QDSnVar,0:PP_N,0:PP_N,1:nSides))
  FluxQDS_Master=0.
  FluxQDS_Slave=0.

  ALLOCATE(QDSMacroValues(1:QDSnVar_macro,0:PP_N,0:PP_N,0:PP_N,nQDSElems))
  QDSMacroValues=0.

  ALLOCATE(GaussHermitWeiAbs(2,2))
  GaussHermitWeiAbs=0.

  ! Weights
  GaussHermitWeiAbs(1,1:2) = 0.8862269254527580136491 ! 0.886227
  GaussHermitWeiAbs(1,:) = GaussHermitWeiAbs(1,:)*SQRT(PI)/SUM(GaussHermitWeiAbs(1,:))
  !Absisc
  GaussHermitWeiAbs(2,1) = -0.7071067811865475244008 !-0.707107
  GaussHermitWeiAbs(2,2) =  0.7071067811865475244008 ! 0.707107
  tempNorm = 0.0
  DO iWeight = 1, 2
    tempNorm = tempNorm + GaussHermitWeiAbs(1,iWeight)*GaussHermitWeiAbs(2,iWeight)**2
  END DO
  GaussHermitWeiAbs(2,:)= GaussHermitWeiAbs(2,:)*SQRT(SQRT(PI)/(2.*tempNorm))

  Temp = 278.687
  Velo =(/0.,0.,0./) !(/1020.882,0.,0./)
  Dens= 2.633459376E25
  Mass = 4.651734101E-26
  QDSSpeciesMass=Mass

  !QDS_Species = GETINT('Particles-QDSSpecies','0')

  QDSSpecDOF=3
  QDSMaxVelo=2*(ABS(Velo(1))+SQRT(2.*BoltzmannConst*Temp/Mass)*ABS(GaussHermitWeiAbs(2,1)))



  !print*,"QDS_dg.f90: ini QDS, do restart?"

  IF(.NOT.DoRestart)THEN
    CALL FillIniQDS()

    !print*,MAXVAL(QDSMacroValues(1,i,j,k,:))
    !print*,MAXVAL(QDSMacroValues(6,i,j,k,:))



!    DO k=0,PP_N
!      DO j=0,PP_N
!        DO i=0,PP_N
!    QDSMacroValues(1,i,j,k,:)=Dens*Mass/10
!    QDSMacroValues(2,i,j,k,:) = QDSMacroValues(1,i,j,k,:)*Velo(1)
!    QDSMacroValues(3,i,j,k,:) = QDSMacroValues(1,i,j,k,:)*Velo(2)
!    QDSMacroValues(4,i,j,k,:) = QDSMacroValues(1,i,j,k,:)*Velo(3)
!    QDSMacroValues(6,i,j,k,:) = Temp
!    END DO; END DO; END DO
!
!    ! fill the cell with no. 88
!    DO k=0,PP_N
!      DO j=0,PP_N
!        DO i=0,PP_N
!    !  i=0; j=0; k=0
!    !  QDSMacroValues(1,i,j,k,1) = Dens*wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,1)*Mass
!      QDSMacroValues(1,i,j,k,88) = Dens*Mass
!      QDSMacroValues(2:4,i,j,k,88) = QDSMacroValues(1,i,j,k,88)*Velo(1:3)
!      QDSMacroValues(6,i,j,k,88) = Temp
!    END DO; END DO; END DO
  END IF

ELSE
  !QDSnVar=0
  QDSnVar_macro=0
END IF
QDSInitDGIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT QDS-DG DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE QDS_InitDG


SUBROUTINE QDSTimeDerivative(t,tStage,tDeriv,doSource,doPrintInfo)
!===================================================================================================================================
! Computes the DG time derivative consisting of Volume Integral and Surface integral for the whole field
! UQDS and UQDSt are allocated
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_QDS_DG_Vars,      ONLY:UQDS,UQDSt,QDSMacroValues,QDSSpeciesMass
USE MOD_QDS_Equation_vars,ONLY:QDSnVar
USE MOD_Vector
USE MOD_QDS_DG_Vars,      ONLY:UQDS,UQDSt,UQDS_master,UQDS_Slave,FluxQDS_Master,FluxQDS_Slave
USE MOD_ProlongToFace,    ONLY:ProlongToFaceQDS
USE MOD_Mesh_Vars,        ONLY:nSides
USE MOD_Interpolation,    ONLY:ApplyJacobianQDS
USE MOD_QDS_SurfInt,      ONLY:SurfIntQDS
USE MOD_QDS_VolInt,       ONLY:VolIntQDS
USE MOD_QDS_FillFlux,     ONLY:FillFluxQDS
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI,              ONLY:StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif /*USE_MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: t,tStage
INTEGER,INTENT(IN)              :: tDeriv
LOGICAL,INTENT(IN)              :: doSource
LOGICAL,INTENT(IN),OPTIONAL     :: doPrintInfo
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(mpiroot.AND.(PRESENT(doPrintInfo)))THEN
  IF(doPrintInfo)THEN
    WRITE(UNIT_stdOut,'(A,ES12.5,A)')' max number density:',MAXVAL(ABS(QDSMacroValues(1,:,:,:,:)))/QDSSpeciesMass,' [1/m^3]'
    IF(MAXVAL(ABS(QDSMacroValues(1,:,:,:,:)))/QDSSpeciesMass.GT.1E50)THEN
      CALL abort(&
      __STAMP__&
      ,'density too high!!! >1e50')
    END IF
  END IF
END IF
! prolong the solution to the face integration points for flux computation
#if USE_MPI
! Prolong to face for MPI sides - send direction
CALL StartReceiveMPIData(QDSnVar,UQDS_Slave,1,nSides,RecRequest_U,SendID=2) ! Receive MINE
CALL ProlongToFaceQDS(UQDS,UQDS_Master,UQDS_Slave,doMPISides=.TRUE.)
CALL StartSendMPIData(QDSnVar,UQDS_Slave,1,nSides,SendRequest_U,SendID=2) ! Send YOUR
#endif /*USE_MPI*/

! Prolong to face for BCSides, InnerSides and MPI sides - receive direction
CALL ProlongToFaceQDS(UQDS,UQDS_Master,UQDS_Slave,doMPISides=.FALSE.)
UQDSt=0.
! compute volume integral contribution and add to ut, first half of all elements
CALL VolIntQDS(UQDSt,dofirstElems=.TRUE.)

#if USE_MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_U,RecRequest_U,SendID=2) !Send YOUR - receive MINE
#endif /*USE_MPI*/

! Initialization of the time derivative
#if USE_MPI
CALL StartReceiveMPIData(QDSnVar,FluxQDS_Slave,1,nSides,RecRequest_Flux,SendID=1) ! Receive MINE
! fill the global surface flux list
CALL FillFluxQDS(t,tDeriv,FluxQDS_Master,FluxQDS_Slave,UQDS_Master,UQDS_Slave,doMPISides=.TRUE.)

CALL StartSendMPIData(QDSnVar,FluxQDS_Slave,1,nSides,SendRequest_Flux,SendID=1) ! Send YOUR
#endif /*USE_MPI*/

! fill the all surface fluxes on this proc
CALL FillFluxQDS(t,tDeriv,FluxQDS_Master,FluxQDS_Slave,UQDS_Master,UQDS_Slave,doMPISides=.FALSE.)
! compute surface integral contribution and add to ut
CALL SurfIntQDS(FluxQDS_Master,FluxQDS_Slave,UQDSt,doMPISides=.FALSE.)

! compute volume integral contribution and add to ut
CALL VolIntQDS(UQDSt,dofirstElems=.FALSE.)

#if USE_MPI
! Complete send / receive
CALL FinishExchangeMPIData(SendRequest_Flux,RecRequest_Flux,SendID=1) !Send MINE -receive YOUR

!FINALIZE Fluxes for MPI Sides
CALL SurfIntQDS(FluxQDS_Master,FluxQDS_Slave,UQDSt,doMPISides=.TRUE.)
#endif

! swap and map to physical space
CALL ApplyJacobianQDS(UQDSt,toPhysical=.TRUE.,toSwap=.TRUE.)

! Add Source Terms
!IF(doSource) CALL CalcSourceQDS(tStage,1.0,Ut)


END SUBROUTINE QDSTimeDerivative


SUBROUTINE FillIniQDS()
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,         ONLY:Elem_xGP
USE MOD_QDS_Equation,      ONLY:QDS_ExactFunc
USE MOD_QDS_DG_Vars,       ONLY:nQDSElems,QDSMacroValues,QDSnVar_macro
USE MOD_QDS_Equation_vars, ONLY:QDSIniExactFunc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iQDSElems
!===================================================================================================================================
DO iQDSElems=1,nQDSElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        CALL QDS_ExactFunc(QDSIniExactFunc,0.,0,Elem_xGP(1:3,i,j,k,iQDSElems),QDSMacroValues(1:QDSnVar_macro,i,j,k,iQDSElems))
      END DO ! i
    END DO ! j
  END DO !k
END DO ! iQDSElems=1,nQDSElems
END SUBROUTINE FillIniQDS



SUBROUTINE QDSReCalculateDGValues()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_QDS_DG_Vars        ,ONLY: QDSSpeciesMass,UQDS,QDSMacroValues,nQDSElems
USE MOD_Mesh_Vars          ,ONLY: sJ
USE MOD_Interpolation_Vars ,ONLY: wGP
USE MOD_QDS_DG_Vars        ,ONLY: QDSnVar_macro
USE MOD_QDS_Equation_vars  ,ONLY: QDSnVar
USE MOD_QDS_Equation       ,ONLY: QDS_Q2U
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem, k, j, i
!===================================================================================================================================
!QDSSpeciesMass=Species(QDS_Species)%MassIC
! Read the maximum number of time steps MaxIter and the end time TEnd from ini file


IF(1.EQ.2)THEN
  DO iElem = 1, nQDSElems
    IF (QDSMacroValues(1,0,0,0,iElem).NE.0.0) then
        print*, 'build parts', iElem
        print*, QDSMacroValues(1,0,0,0,iElem)/(QDSSpeciesMass*wGP(0)*wGP(0)*wGP(0)/sJ(0,0,0,iElem)), &
        QDSMacroValues(2:4,0,0,0,iElem)/QDSMacroValues(1,0,0,0,iElem), QDSMacroValues(6,0,0,0,iElem)
    end if
    IF (QDSMacroValues(1,1,1,1,iElem).NE.0.0) then
        print*, 'zwei', iElem
        print*, QDSMacroValues(1,1,1,1,iElem)/(QDSSpeciesMass*wGP(1)*wGP(1)*wGP(1)/sJ(1,1,1,iElem)), &
        QDSMacroValues(2:4,1,1,1,iElem)/QDSMacroValues(1,1,1,1,iElem), QDSMacroValues(6,1,1,1,iElem)
    end if
  END DO
END IF


DO iElem = 1, nQDSElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        !L = 0
        IF (QDSMacroValues(1,i,j,k,iElem).GT.0.0) THEN
          IF (QDSMacroValues(6,i,j,k,iElem).LT.0.0) QDSMacroValues(6,i,j,k,iElem) = 0.0
            CALL QDS_Q2U(QDSMacroValues(1:QDSnVar_macro,i,j,k,iElem), UQDS(1:QDSnVar,i,j,k,iElem))

          !DO iPart1=1,2; DO iPart2=1,2; DO iPart3=1,2
          !  UQDS(1+L,i,j,k,iElem) =     QDSMacroValues(1,i,j,k,iElem)*&
          !                           GaussHermitWeiAbs(1,iPart1)     *&
          !                           GaussHermitWeiAbs(1,iPart2)     *&
          !                           GaussHermitWeiAbs(1,iPart3)/(PI*SQRT(PI))

          !  UQDS(2+L,i,j,k,iElem) =               UQDS(1+L,i,j,k,iElem) &
          !                           * (QDSMacroValues(2  ,i,j,k,iElem) /&
          !                              QDSMacroValues(1  ,i,j,k,iElem) &
          !       + SQRT(2.*BoltzmannConst*QDSMacroValues(6,i,j,k,iElem)/QDSSpeciesMass)*GaussHermitWeiAbs(2,iPart1))

          !  UQDS(3+L,i,j,k,iElem) =               UQDS(1+L,i,j,k,iElem) &
          !                           * (QDSMacroValues(3  ,i,j,k,iElem) /&
          !                              QDSMacroValues(1  ,i,j,k,iElem) &
          !       + SQRT(2.*BoltzmannConst*QDSMacroValues(6,i,j,k,iElem)/QDSSpeciesMass)*GaussHermitWeiAbs(2,iPart2))

          !  UQDS(4+L,i,j,k,iElem) =               UQDS(1+L,i,j,k,iElem) &
          !                           * (QDSMacroValues(4  ,i,j,k,iElem) /&
          !                              QDSMacroValues(1  ,i,j,k,iElem) &
          !       + SQRT(2.*BoltzmannConst*QDSMacroValues(6,i,j,k,iElem)/QDSSpeciesMass)*GaussHermitWeiAbs(2,iPart3))

          !  UQDS(5+L,i,j,k,iElem) =(QDSSpecDOF-3.)*BoltzmannConst*QDSMacroValues(6,i,j,k,iElem)/(QDSSpeciesMass*2.)

          !  L = L + 5
          !END DO; END DO; END DO
        ELSE
          UQDS(:,i,j,k,iElem) = 0.0
        END IF
      END DO ! i
    END DO ! j
  END DO
END DO
END SUBROUTINE QDSReCalculateDGValues


SUBROUTINE QDSCalculateMacroValues()
!===================================================================================================================================
! Get the constant advection velocity vector from the ini file
!===================================================================================================================================
! MODULES
USE MOD_QDS_DG_Vars
USE MOD_Globals_Vars,       ONLY:BoltzmannConst
!USE MOD_QDS_DG_Vars,           ONLY : QDS_Species
USE MOD_QDS_DG_Vars,           ONLY : QDSSpeciesMass
USE MOD_PreProc

!USE MOD_Mesh_Vars,          ONLY : sJ
!USE MOD_Interpolation_Vars, ONLY : wGP
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem, k, j, i, iPart,L
!REAL :: Temp, Velo(3), Dens, Mass
!===================================================================================================================================
! Read the maximum number of time steps MaxIter and the end time TEnd from ini file
DO iElem = 1, nQDSElems
  DO k=0,PP_N
    DO j=0,PP_N
      DO i=0,PP_N
        QDSMacroValues(:,i,j,k,iElem) = 0.0
        DO iPart=0,7
          L=iPart*5
          IF(UQDS(1+L,i,j,k,iElem).GT.0.0)THEN
            QDSMacroValues(1,i,j,k,iElem) = QDSMacroValues(1,i,j,k,iElem) + UQDS(1+L,i,j,k,iElem)
            QDSMacroValues(2,i,j,k,iElem) = QDSMacroValues(2,i,j,k,iElem) + UQDS(2+L,i,j,k,iElem)
            QDSMacroValues(3,i,j,k,iElem) = QDSMacroValues(3,i,j,k,iElem) + UQDS(3+L,i,j,k,iElem)
            QDSMacroValues(4,i,j,k,iElem) = QDSMacroValues(4,i,j,k,iElem) + UQDS(4+L,i,j,k,iElem)
            QDSMacroValues(5,i,j,k,iElem) = QDSMacroValues(5,i,j,k,iElem) + UQDS(5+L,i,j,k,iElem) &
                                                                          + 0.5*(UQDS(2+L,i,j,k,iElem)**2+&
                                                                                 UQDS(3+L,i,j,k,iElem)**2+&
                                                                                 UQDS(4+L,i,j,k,iElem)**2)&
                                                                                /UQDS(1+L,i,j,k,iElem)
          END IF
        END DO
        IF (QDSMacroValues(1,i,j,k,iElem).GT.0.0) THEN
          QDSMacroValues(6,i,j,k,iElem) =  (2.*QDSMacroValues(5,i,j,k,iElem) &
                                             -(QDSMacroValues(2,i,j,k,iElem)**2+&
                                               QDSMacroValues(3,i,j,k,iElem)**2+&
                                               QDSMacroValues(4,i,j,k,iElem)**2) &
            /QDSMacroValues(1,i,j,k,iElem)) / (QDSMacroValues(1,i,j,k,iElem)*3.) *QDSSpeciesMass /BoltzmannConst
        ELSE
          QDSMacroValues(1:6,i,j,k,iElem) = 0.0
        END IF
      END DO ! i
    END DO ! j
  END DO
END DO
END SUBROUTINE QDSCalculateMacroValues


SUBROUTINE QDS_FinalizeDG()
!===================================================================================================================================
!>
!===================================================================================================================================
! MODULES
USE MOD_QDS_DG_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
QDSInitDGIsDone = .FALSE.
SDEALLOCATE(UQDS)
SDEALLOCATE(UQDSt)
SDEALLOCATE(UQDS_Master)
SDEALLOCATE(UQDS_Slave)
SDEALLOCATE(FluxQDS_Master)
SDEALLOCATE(FluxQDS_Slave)
SDEALLOCATE(GaussHermitWeiAbs)
END SUBROUTINE QDS_FinalizeDG


END MODULE MOD_QDS_DG
#endif /*USE_QDS_DG*/
