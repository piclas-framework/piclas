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

MODULE MOD_DSMC_SteadyState
!===================================================================================================================================
! Module for DSMC
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

INTERFACE DSMC_SteadyStateInit
  MODULE PROCEDURE DSMC_SteadyStateInit
END INTERFACE

INTERFACE QCrit_evaluation
  MODULE PROCEDURE QCrit_evaluation
END INTERFACE

INTERFACE SteadyStateDetection_main
  MODULE PROCEDURE SteadyStateDetection_main
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC :: DSMC_SteadyStateInit, QCrit_evaluation, SteadyStateDetection_main
PUBLIC :: SteadyStateDetection_Sampling, EntropyCalculation, SteadyStateDetection_Algorithm, KernelDensityEstimation
!===================================================================================================================================

CONTAINS

SUBROUTINE DSMC_SteadyStateInit()
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
  USE MOD_Globals
  USE MOD_Mesh_Vars,             ONLY : nElems, nBCSides
  USE MOD_ReadInTools
  USE MOD_DSMC_Vars
  USE MOD_PARTICLE_Vars,         ONLY: nSpecies
! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

  IF (UseQCrit) THEN
    SWRITE(*,*)'Steady State Detection with Q-Criterion'
    QCritTestStep = GETINT('Particles-DSMC-QCritTestStep','1000')
    QCritEpsilon  = GETREAL('Particles-DSMC-QCritEpsilon','0.1')
    ALLOCATE(QCritCounter(nBCSides,2))
    QCritCounter(1:nBCSides,1:2) = 0
    ALLOCATE(QLocal(nElems))
    QLocal(1:nElems) = -99
    QCritLastTest = 0
  ELSEIF(UseSSD) THEN
    SWRITE(*,*)'Steady State Detection with SSD-Algorithm'
    iSamplingIters = 0
    nSamplingIters = GETINT('Particles-DSMC-SSD-SampIters','100')
    HistTime = 0
    nTime = GETINT('Particles-DSMC-SSD-HistSize','30')
    ! Set Parameters for the Von Neumann Ratio - Test
    ALLOCATE(RValue(nElems))
    RValue(1:nElems) = 0
    SELECT CASE (nTime)
     CASE (10)
      !Critical(1) = 0.5312   ! 10% Significance
      !Critical(2) = 1.4689   ! ""
      Critical(1) = 0.3759   ! 2% Significance
      Critical(2) = 1.6241   ! ""
     CASE (20)
      !Critical(1) = 0.6498
      !Critical(2) = 1.3502
      Critical(1) = 0.5200
      Critical(2) = 1.4797
     CASE (30)
      !Critical(1) = 0.7092
      !Critical(2) = 1.2909
      Critical(1) = 0.5976
      Critical(2) = 1.4025
     CASE (40)
      !Critical(1) = 0.7461
      !Critical(2) = 1.2540
      Critical(1) = 0.6467
      Critical(2) = 1.3533
     CASE (50)
      !Critical(1) = 0.7719
      !Critical(2) = 1.2282
      Critical(1) = 0.6815
      Critical(2) = 1.3185
     CASE (60)
      !Critical(1) = 0.7907
      !Critical(2) = 1.2093
      Critical(1) = 0.7072
      Critical(2) = 1.2928
     CASE DEFAULT
      Critical(1) = 0.5200
      Critical(2) = 1.4797
      SWRITE(*,*)'Time Step not available!',nTime
    END SELECT
    ! Set Parameters for the Euclidian Distance Criterion
    ALLOCATE(ED_Delta(nElems))
    ED_Delta(1:nElems) = -2
    SELECT CASE (nTime)
     CASE (5)
      !Epsilon1 = 0.95110  ! 5% Significance
      !Epsilon2 = 1.9994  ! ""
      Epsilon1 = 1.15300  ! 2% Significance
      Epsilon2 = 1.9996    ! ""
     CASE (8)
      !Epsilon1 = 0.86910
      !Epsilon2 = 1.9991
      Epsilon1 = 1.04100
      Epsilon2 = 1.9992
     CASE (10)
      !Epsilon1 = 0.850645
      !Epsilon2 = 2.0032
      Epsilon1 = 0.97980
      Epsilon2 = 1.9991
     CASE (15)
      !Epsilon1 = 0.79284
      !Epsilon2 = 2.0032
      Epsilon1 = 0.86091
      Epsilon2 = 1.9906
     CASE (20)
      !Epsilon1 = 0.75450
      !Epsilon2 = 1.9991
      Epsilon1 = 0.81049
      Epsilon2 = 1.9906
     CASE (25)
      !Epsilon1 = 0.73164
      !Epsilon2 = 2.0032
      Epsilon1 = 0.77867
      Epsilon2 = 1.9907
     CASE (30)
      !Epsilon1 = 0.68313
      !Epsilon2 = 1.9906
      Epsilon1 = 0.74786
      Epsilon2 = 1.9906
     CASE (35)
      Epsilon1 = 0.68270
      Epsilon2 = 1.9990
      !Epsilon1 = 0.76210
      !Epsilon2 = 1.9990
     CASE DEFAULT
      Epsilon1 = 0.81049
      Epsilon2 = 1.9906
      SWRITE(*,*)'Time Step not available!',nTime
    END SELECT
    ! Set Parameters for the Student-t-Test
    ALLOCATE(Stud_Indicator(nElems))
    Stud_Indicator(1:nElems) = 0
    SELECT CASE (nTime)
     CASE (5)
      !StudCrit = 4.011    ! 5% Significance
      StudCrit = 5.021     ! 2% Significance
     CASE (10)
      !StudCrit = 3.568
      StudCrit = 4.138
     CASE (15)
      !StudCrit = 3.472
      StudCrit = 3.927
     CASE (20)
      !StudCrit = 3.445
      StudCrit = 3.845
     CASE (25)
      !StudCrit = 3.440
      StudCrit = 3.809
     CASE (30)
      !StudCrit = 3.445
      StudCrit = 3.792
     CASE (40)
      !StudCrit = 3.464
      StudCrit = 3.785
     CASE DEFAULT
      StudCrit = 3.845
      SWRITE(*,*)'Time Step not available!',nTime
    END SELECT
    ! Set Parameters for the Polynomial Interpolation Test
    ALLOCATE(ConvCoeff(nTime))
    ConvCoeff(1:nTime) = 0
    ALLOCATE(PIT_Drift(nElems))
    PIT_Drift(1:nElems) = 2
    SELECT CASE (nTime)
     CASE (5)
      !PITCrit = 1.862  !C = 0.95, Cubic, NV  ! 5% Significance
      !PITCrit = 2.210  !C = 0.98, Cubic, NV  ! 2% Significance
      !PITCrit = 2.447  !C = 0.99, Cubic, NV  ! 1% Significance
      !PITCrit = 2.443  !C = 0.95, Cubic, StudT  ! 5% Significance
      PITCrit = 3.197  !C = 0.98, Cubic, StudT  ! 2% Significance
      !PITCrit = 3.831  !C = 0.99, Cubic, StudT  ! 1% Significance
      !Cubic Interpolation
      ConvCoeff(1) = 1.0/12
      ConvCoeff(2) = -8.0/12
      ConvCoeff(3) = 0
      ConvCoeff(4) = - ConvCoeff(2)
      ConvCoeff(5) = - ConvCoeff(1)
     CASE (7)
      !PITCrit = 1.004
      !PITCrit = 1.192
      !PITCrit = 1.320
      !PITCrit = 1.212
      PITCrit = 1.536
      !PITCrit = 1.793
      !Cubic Interpolation
      ConvCoeff(1) = 22.0/252
      ConvCoeff(2) = -67.0/252
      ConvCoeff(3) = -58.0/252
      ConvCoeff(4) = 0
      ConvCoeff(5) = - ConvCoeff(3)
      ConvCoeff(6) = - ConvCoeff(2)
      ConvCoeff(7) = - ConvCoeff(1)
     CASE (9)
      !PITCrit = 0.6628
      !PITCrit = 0.7865
      !PITCrit = 0.8710
      !PITCrit = 0.7649
      PITCrit = 0.9539
      !PITCrit = 1.099
      !Cubic Interpolation
      ConvCoeff(1) = 86.0/1188
      ConvCoeff(2) = -142.0/1188
      ConvCoeff(3) = -193.0/1188
      ConvCoeff(4) = -126.0/1188
      ConvCoeff(5) = 0
      ConvCoeff(6) = - ConvCoeff(4)
      ConvCoeff(7) = - ConvCoeff(3)
      ConvCoeff(8) = - ConvCoeff(2)
      ConvCoeff(9) = - ConvCoeff(1)
     CASE (11)
      !PITCrit = 0.4816
      !PITCrit = 0.5715
      !PITCrit = 0.6330
      !PITCrit = 0.5408
      PITCrit = 0.6679
      !PITCrit = 0.7632
      !Cubic Interpolation
      ConvCoeff(1) = 300.0/5148
      ConvCoeff(2) = -294.0/5148
      ConvCoeff(3) = -532.0/5148
      ConvCoeff(4) = -503.0/5148
      ConvCoeff(5) = -296.0/5148
      ConvCoeff(6) = 0
      ConvCoeff(7) = - ConvCoeff(5)
      ConvCoeff(8) = - ConvCoeff(4)
      ConvCoeff(9) = - ConvCoeff(3)
      ConvCoeff(10) = - ConvCoeff(2)
      ConvCoeff(11) = - ConvCoeff(1)
     CASE (13)
      !PITCrit = 0.3711
      !PITCrit = 0.4404
      !PITCrit = 0.4877
      !PITCrit = 0.4090
      PITCrit = 0.5017
      !PITCrit = 0.5703
      !Cubic Interpolation
      ConvCoeff(1) = 1133.0/24024
      ConvCoeff(2) = -660.0/24024
      ConvCoeff(3) = -1578.0/24024
      ConvCoeff(4) = -1796.0/24024
      ConvCoeff(5) = -1489.0/24024
      ConvCoeff(6) = -832.0/24024
      ConvCoeff(7) = 0
      ConvCoeff(8) = - ConvCoeff(6)
      ConvCoeff(9) = - ConvCoeff(5)
      ConvCoeff(10) = - ConvCoeff(4)
      ConvCoeff(11) = - ConvCoeff(3)
      ConvCoeff(12) = - ConvCoeff(2)
      ConvCoeff(13) = - ConvCoeff(1)
     CASE (15)
      !PITCrit = 0.2975
      !PITCrit = 0.3531
      !PITCrit = 0.3911
      !PITCrit = 0.3235
      PITCrit = 0.3950
      !PITCrit = 0.4474
      !Cubic Interpolation
      ConvCoeff(1) = 12922.0/334152
      ConvCoeff(2) = -4121.0/334152
      ConvCoeff(3) = -14150.0/334152
      ConvCoeff(4) = -18334.0/334152
      ConvCoeff(5) = -17842.0/334152
      ConvCoeff(6) = -13843.0/334152
      ConvCoeff(7) = -7506.0/334152
      ConvCoeff(8) = 0
      ConvCoeff(9) = - ConvCoeff(7)
      ConvCoeff(10) = - ConvCoeff(6)
      ConvCoeff(11) = - ConvCoeff(5)
      ConvCoeff(12) = - ConvCoeff(4)
      ConvCoeff(13) = - ConvCoeff(3)
      ConvCoeff(14) = - ConvCoeff(2)
      ConvCoeff(15) = - ConvCoeff(1)
     CASE (17)
      !PITCrit = 0.2456
      !PITCrit = 0.2915
      !PITCrit = 0.3228
      !PITCrit = 0.2644
      PITCrit = 0.3217
      !PITCrit = 0.3631
      !Cubic Interpolation
      ConvCoeff(1) = 748.0/23256
      ConvCoeff(2) = -98.0/23256
      ConvCoeff(3) = -643.0/23256
      ConvCoeff(4) = -930.0/23256
      ConvCoeff(5) = -1002.0/23256
      ConvCoeff(6) = -902.0/23256
      ConvCoeff(7) = -673.0/23256
      ConvCoeff(8) = -358.0/23256
      ConvCoeff(9) = 0
      ConvCoeff(10) = - ConvCoeff(8)
      ConvCoeff(11) = - ConvCoeff(7)
      ConvCoeff(12) = - ConvCoeff(6)
      ConvCoeff(13) = - ConvCoeff(5)
      ConvCoeff(14) = - ConvCoeff(4)
      ConvCoeff(15) = - ConvCoeff(3)
      ConvCoeff(16) = - ConvCoeff(2)
      ConvCoeff(17) = - ConvCoeff(1)
     CASE (19)
      !PITCrit = 0.2073
      !PITCrit = 0.2460
      !PITCrit = 0.2724
      !PITCrit = 0.2213
      PITCrit = 0.2685
      !PITCrit = 0.3026
      !Cubic Interpolation
      ConvCoeff(1) = 6936.0/255816
      ConvCoeff(2) = 68.0/255816
      ConvCoeff(3) = -4648.0/255816
      ConvCoeff(4) = -7481.0/255816
      ConvCoeff(5) = -8700.0/255816
      ConvCoeff(6) = -8574.0/255816
      ConvCoeff(7) = -7372.0/255816
      ConvCoeff(8) = -5363.0/255816
      ConvCoeff(9) = -2816.0/255816
      ConvCoeff(10) = 0
      ConvCoeff(11) = - ConvCoeff(9)
      ConvCoeff(12) = - ConvCoeff(8)
      ConvCoeff(13) = - ConvCoeff(7)
      ConvCoeff(14) = - ConvCoeff(6)
      ConvCoeff(15) = - ConvCoeff(5)
      ConvCoeff(16) = - ConvCoeff(4)
      ConvCoeff(17) = - ConvCoeff(3)
      ConvCoeff(18) = - ConvCoeff(2)
      ConvCoeff(19) = - ConvCoeff(1)
     CASE (21)
      !PITCrit = 0.1780
      !PITCrit = 0.2112
      !PITCrit = 0.2340
      !PITCrit = 0.1889
      PITCrit = 0.2287
      !PITCrit = 0.2571
      !Cubic Interpolation
      ConvCoeff(1) = 84075.0/3634092
      ConvCoeff(2) = 10032.0/3634092
      ConvCoeff(3) = -43284.0/3634092
      ConvCoeff(4) = -78176.0/3634092
      ConvCoeff(5) = -96947.0/3634092
      ConvCoeff(6) = -101900.0/3634092
      ConvCoeff(7) = -95338.0/3634092
      ConvCoeff(8) = -79564.0/3634092
      ConvCoeff(9) = -56881.0/3634092
      ConvCoeff(10) = -29592.0/3634092
      ConvCoeff(11) = 0
      ConvCoeff(12) = - ConvCoeff(10)
      ConvCoeff(13) = - ConvCoeff(9)
      ConvCoeff(14) = - ConvCoeff(8)
      ConvCoeff(15) = - ConvCoeff(7)
      ConvCoeff(16) = - ConvCoeff(6)
      ConvCoeff(17) = - ConvCoeff(5)
      ConvCoeff(18) = - ConvCoeff(4)
      ConvCoeff(19) = - ConvCoeff(3)
      ConvCoeff(20) = - ConvCoeff(2)
      ConvCoeff(21) = - ConvCoeff(1)
     CASE (23)
      !PITCrit = 0.1551
      !PITCrit = 0.1840
      !PITCrit = 0.2038
      !PITCrit = 0.1637
      PITCrit = 0.1978
      !PITCrit = 0.2221
      !Cubic Interpolation
      ConvCoeff(1) = 3938.0/197340
      ConvCoeff(2) = 815.0/197340
      ConvCoeff(3) = -1518.0/197340
      ConvCoeff(4) = -3140.0/197340
      ConvCoeff(5) = -4130.0/197340
      ConvCoeff(6) = -4567.0/197340
      ConvCoeff(7) = -4530.0/197340
      ConvCoeff(8) = -4098.0/197340
      ConvCoeff(9) = -3350.0/197340
      ConvCoeff(10) = -2365.0/197340
      ConvCoeff(11) = -1222.0/197340
      ConvCoeff(12) = 0
      ConvCoeff(13) = - ConvCoeff(11)
      ConvCoeff(14) = - ConvCoeff(10)
      ConvCoeff(15) = - ConvCoeff(9)
      ConvCoeff(16) = - ConvCoeff(8)
      ConvCoeff(17) = - ConvCoeff(7)
      ConvCoeff(18) = - ConvCoeff(6)
      ConvCoeff(19) = - ConvCoeff(5)
      ConvCoeff(20) = - ConvCoeff(4)
      ConvCoeff(21) = - ConvCoeff(3)
      ConvCoeff(22) = - ConvCoeff(2)
      ConvCoeff(23) = - ConvCoeff(1)
     CASE (25)
      !PITCrit = 0.1367
      !PITCrit = 0.1622
      !PITCrit = 0.1796
      !PITCrit = 0.1436
      PITCrit = 0.1733
      !PITCrit = 0.1943
      !Cubic Interpolation
      ConvCoeff(1) = 30866.0/1776060
      ConvCoeff(2) = 8602.0/1776060
      ConvCoeff(3) = -8525.0/1776060
      ConvCoeff(4) = -20982.0/1776060
      ConvCoeff(5) = -29236.0/1776060
      ConvCoeff(6) = -33754.0/1776060
      ConvCoeff(7) = -35003.0/1776060
      ConvCoeff(8) = -33450.0/1776060
      ConvCoeff(9) = -29562.0/1776060
      ConvCoeff(10) = -23806.0/1776060
      ConvCoeff(11) = -16649.0/1776060
      ConvCoeff(12) = -8558.0/1776060
      ConvCoeff(13) = 0
      ConvCoeff(14) = - ConvCoeff(12)
      ConvCoeff(15) = - ConvCoeff(11)
      ConvCoeff(16) = - ConvCoeff(10)
      ConvCoeff(17) = - ConvCoeff(9)
      ConvCoeff(18) = - ConvCoeff(8)
      ConvCoeff(19) = - ConvCoeff(7)
      ConvCoeff(20) = - ConvCoeff(6)
      ConvCoeff(21) = - ConvCoeff(5)
      ConvCoeff(22) = - ConvCoeff(4)
      ConvCoeff(23) = - ConvCoeff(3)
      ConvCoeff(24) = - ConvCoeff(2)
      ConvCoeff(25) = - ConvCoeff(1)
     CASE DEFAULT
      SWRITE(*,*)'Time Step not available!',nTime
    END SELECT
  ! Set Parameters for the Mann - Kendall - Test
    ALLOCATE(MK_Trend(nElems))
    MK_Trend(1:nElems) = 10

    ALLOCATE(Sampler(nElems,nSpecies))
    Sampler(1:nElems,1:nSpecies)%Energy(1) = 0.0
    Sampler(1:nElems,1:nSpecies)%Energy(2) = 0.0
    Sampler(1:nElems,1:nSpecies)%Energy(3) = 0.0
    Sampler(1:nElems,1:nSpecies)%Velocity(1) = 0.0
    Sampler(1:nElems,1:nSpecies)%Velocity(2) = 0.0
    Sampler(1:nElems,1:nSpecies)%Velocity(3) = 0.0
    Sampler(1:nElems,1:nSpecies)%PartNum = 0.0
    Sampler(1:nElems,1:nSpecies)%ERot = 0.0
    Sampler(1:nElems,1:nSpecies)%EVib = 0.0
    Sampler(1:nElems,1:nSpecies)%EElec = 0.0
    ALLOCATE(History(nElems,nSpecies,nTime))
    History(1:nElems,1:nSpecies,1:nTime)%Energy(1) = 0.0
    History(1:nElems,1:nSpecies,1:nTime)%Energy(2) = 0.0
    History(1:nElems,1:nSpecies,1:nTime)%Energy(3) = 0.0
    History(1:nElems,1:nSpecies,1:nTime)%Velocity(1) = 0.0
    History(1:nElems,1:nSpecies,1:nTime)%Velocity(2) = 0.0
    History(1:nElems,1:nSpecies,1:nTime)%Velocity(3) = 0.0
    History(1:nElems,1:nSpecies,1:nTime)%PartNum = 0.0
    History(1:nElems,1:nSpecies,1:nTime)%ERot = 0.0
    History(1:nElems,1:nSpecies,1:nTime)%EVib = 0.0
    History(1:nElems,1:nSpecies,1:nTime)%EElec = 0.0
    ALLOCATE(CheckHistory(nElems,nTime))
    CheckHistory(1:nElems,1:nTime) = 0.0
    ALLOCATE(SteadyIdent(nElems,nSpecies,10))
    SteadyIdent(1:nElems,1:nSpecies,1:10) = 0
    ALLOCATE(SteadyIdentGlobal(nSpecies,10))
    SteadyIdentGlobal(1:nSpecies,1:10) = 0
  ENDIF

END SUBROUTINE DSMC_SteadyStateInit


SUBROUTINE  QCrit_evaluation()

  USE MOD_DSMC_Vars,             ONLY : SamplingActive, QCritCounter, QCritEpsilon
  USE MOD_Mesh_Vars,             ONLY : nBCSides, BC
  USE MOD_Particle_Boundary_Vars,ONLY : PartBound
  USE MOD_TimeDisc_Vars,         ONLY : iter
  USE MOD_DSMC_Vars,             ONLY : QLocal
#ifdef MPI
  USE mpi
  USE MOD_Globals,               ONLY : MPIRoot
#endif
!--------------------------------------------------------------------------------------------------!
! Q-Criterion Evaluation Routine (Boyd,Burt)
!  Count the number of particles that interact with surfaces or domain boundaries per cell
!  over 'QCritTestStep' iterations (particle_boundary_treatment).
!  Analyze the difference between two successive intervals statistically for steady state conditions.
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  REAL              :: Qfactor1                   ! Factors for Normation in Q-Value Q = max * 1/N
  REAL              :: Qfactor2                   !  N = sqrt(Qfactor1 + Qfactor2 * log(nSides))
  REAL              :: QValue
  REAL              :: maxValue, SideValue
  INTEGER           :: iSide, nSides

  INTEGER :: rank ! DEBUG
#ifdef MPI
  REAL              :: maxValue_global
  INTEGER           :: nSides_global
  INTEGER           :: IERROR
#endif
!--------------------------------------------------------------------------------------------------!

  Qfactor1 = -1.23 ! Boyd, Burt
  Qfactor2 = 1.85  ! Boyd, Burt

  ! DEBUG
#ifdef MPI
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,IERROR)
#else
  rank = 0
#endif
  ! DEBUG

  maxValue = 0.0
  nSides = 0
  DO iSide=1,nBCSides
    IF(((QCritCounter(iSide,1)+QCritCounter(iSide,2)).GT.0).AND.PartBound%UseForQCrit(BC(iSide))) THEN
      ! Evaluate Criterion for side iSide with successive counts:
      SideValue = abs(REAL(QCritCounter(iSide,1)-QCritCounter(iSide,2)))/sqrt(REAL(QCritCounter(iSide,1)+QCritCounter(iSide,2)))
      QLocal(iSide) = SideValue
      IF(SideValue.GE.maxValue) maxValue = SideValue
      ! Save latest Count, Reset acutal Counter
      QCritCounter(iSide,2) = QCritCounter(iSide,1)
      QCritCounter(iSide,1) = 0
      nSides = nSides + 1
    ENDIF
  ENDDO

#ifdef MPI
      CALL MPI_ALLREDUCE(nSides,nSides_global,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
      nSides = nSides_global
#endif

  QLocal(1:nBCSides) = QLocal(1:nBCSides) / sqrt(Qfactor1+Qfactor2*log(REAL(nSides)))

  IF(nSides.GE.20) THEN
#ifdef MPI
    CALL MPI_ALLREDUCE(maxValue,maxValue_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    maxValue = maxValue_global
#endif
    QValue = maxValue / sqrt(Qfactor1+Qfactor2*log(REAL(nSides)))
  ELSE
    SWRITE(*,*)'ERROR: Not enough Sides for Q-Criterion'
    QValue = 0.0
  ENDIF

! DEBUG
  IF(rank.EQ.0) THEN
      OPEN(UNIT=900,FILE='QCrit.txt',POSITION='APPEND')
      WRITE(900,*)'QCrit Test: Iteration Q-Value',iter+1,QValue
      CLOSE(900)
  ENDIF
! DEBUG

  IF(QValue.LT.(1.0+QCritEpsilon)) SamplingActive = .TRUE.

END SUBROUTINE QCrit_evaluation


SUBROUTINE SteadyStateDetection_main()

  USE MOD_DSMC_Vars,             ONLY : DSMC, CollisMode, SpecDSMC
  USE MOD_Mesh_Vars,             ONLY : nElems
  USE MOD_Particle_Vars,         ONLY : nSpecies, Species
  USE MOD_Globals_Vars,          ONLY : BoltzmannConst
  USE MOD_DSMC_Vars,             ONLY : iSamplingIters, nSamplingIters, HistTime, nTime
  USE MOD_DSMC_Vars,             ONLY : Sampler, History, CheckHistory, SamplingActive, SteadyIdentGlobal
#ifdef MPI
  USE mpi
#endif

!--------------------------------------------------------------------------------------------------------!
! Steady State Detection - Algorithm
!  1. Sampling of local flow properties over 'nSamplingIters' for variance reduction
!  2. Analyze the last 'nTime' averaged values for steady state conditions using a statistical test
!   - all properties+species combinations are treated separately
!   - testing each cell separately, then derive the global state
!--------------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                         !
!--------------------------------------------------------------------------------------------------------!
! argument list declaration                                                                              !
! Local variable declaration
  INTEGER            :: iElem,iSpecies,iTime,iVal
  LOGICAL            :: CheckVelo(3)   ! Flagg whether flow property is to check: velocity (x,y,z)
  LOGICAL            :: CheckEnerg(3)  !  trans. temperature (x,y,z)
  LOGICAL            :: CheckDens  !  density
  LOGICAL            :: CheckIntEn(3)  !  rotational energy, vibrational energy, electronic state
  INTEGER            :: Check(10)  ! Number of Species identified as 'steady state' for every property of the above,
          !  i.e. if Check(i)=nSpecies, then property i is in steady state
  INTEGER            :: TotCheck        ! Number of flow properties that are in steady state
!--------------------------------------------------------------------------------------------------------!

! Sampling of macroscopic properties (over nSamplingIters) for reduction of variance and autocorrelation effects (Step 1)
  CALL SteadyStateDetection_Sampling()

  ! initialize all flow properties to check with false
  CheckVelo = .FALSE.
  CheckEnerg= .FALSE.
  CheckIntEn= .FALSE.
  CheckDens = .FALSE.
  ! Set the flow properties to check
  CheckVelo(1) = .TRUE.  ! x - velocity
  !CheckVelo(2) = .TRUE.  ! y - velocity
  !CheckVelo(3) = .TRUE.  ! z - velocity
  !CheckEnerg(1) = .TRUE. ! x - trans. temperature
  !CheckEnerg(2) = .TRUE. ! y - trans. temperature
  !CheckEnerg(3) = .TRUE. ! z - trans. temperature
  !CheckIntEn(1) = .TRUE. ! rotational energy
  !CheckIntEn(2) = .TRUE. ! vibrational energy
  !CheckIntEn(3) = .TRUE. ! electronic state
  !CheckDens = .TRUE.   ! density

  iSamplingIters = iSamplingIters + 1

  IF(iSamplingIters.EQ.nSamplingIters) THEN
! Step 2: Write Sampled Values into History (memorizing last nTime sampled values only)
    IF(HistTime.LT.nTime) THEN
      HistTime = HistTime + 1
    ELSE
      DO iSpecies = 1, nSpecies
        DO iElem = 1, nElems
          DO iTime = 1, nTime-1
            History(iElem,iSpecies,iTime)%Energy(1:3)     = History(iElem,iSpecies,iTime+1)%Energy(1:3)
            History(iElem,iSpecies,iTime)%Velocity(1:3)   = History(iElem,iSpecies,iTime+1)%Velocity(1:3)
            History(iElem,iSpecies,iTime)%PartNum         = History(iElem,iSpecies,iTime+1)%PartNum
            IF((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
              IF ((SpecDSMC(iSpecies)%InterID.EQ.2).OR.(SpecDSMC(iSpecies)%InterID.EQ.20)) THEN
                History(iElem,iSpecies,iTime)%ERot          = History(iElem,iSpecies,iTime+1)%ERot
                History(iElem,iSpecies,iTime)%EVib          = History(iElem,iSpecies,iTime+1)%EVib
              END IF
            END IF
            IF(DSMC%ElectronicModel) THEN
              History(iElem,iSpecies,iTime)%EElec         = History(iElem,iSpecies,iTime+1)%EElec
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    DO iSpecies = 1, nSpecies
      DO iElem = 1, nElems
        IF(Sampler(iElem,iSpecies)%PartNum.GT.0) THEN
          History(iElem,iSpecies,HistTime)%Velocity(1:3) = Sampler(iElem,iSpecies)%Velocity(1:3) &
                                                         / Sampler(iElem,iSpecies)%PartNum
          History(iElem,iSpecies,HistTime)%Energy(1:3)   = Sampler(iElem,iSpecies)%Energy(1:3) / Sampler(iElem,iSpecies)%PartNum &
                                                         - 0.5 * History(iElem,iSpecies,HistTime)%Velocity(1:3)**2
          ! Show as translational temperature
          History(iElem,iSpecies,HistTime)%Energy(1:3)   = 2 * History(iElem,iSpecies,HistTime)%Energy(1:3) &
                                                         * Species(iSpecies)%MassIC &
                                                         / BoltzmannConst
          History(iElem,iSpecies,HistTime)%PartNum       = Sampler(iElem,iSpecies)%PartNum
          IF((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
            IF((SpecDSMC(iSpecies)%InterID.EQ.2).OR.(SpecDSMC(iSpecies)%InterID.EQ.20)) THEN
              History(iElem,iSpecies,HistTime)%ERot        = Sampler(iElem,iSpecies)%ERot  / Sampler(iElem,iSpecies)%PartNum
              History(iElem,iSpecies,HistTime)%EVib        = Sampler(iElem,iSpecies)%EVib  / Sampler(iElem,iSpecies)%PartNum
            END IF
          END IF
          IF(DSMC%ElectronicModel) THEN
            History(iElem,iSpecies,HistTime)%EElec       = Sampler(iElem,iSpecies)%EElec / Sampler(iElem,iSpecies)%PartNum
          ENDIF
        ELSE
          History(iElem,iSpecies,HistTime)%Velocity(1:3) = -1.0
          History(iElem,iSpecies,HistTime)%Energy(1:3)   = -1.0
          History(iElem,iSpecies,HistTime)%PartNum       = -1.0
          IF((CollisMode.EQ.2).OR.(CollisMode.EQ.3))THEN
            IF ((SpecDSMC(iSpecies)%InterID.EQ.2).OR.(SpecDSMC(iSpecies)%InterID.EQ.20)) THEN
              History(iElem,iSpecies,HistTime)%ERot        = -1.0
              History(iElem,iSpecies,HistTime)%EVib        = -1.0
            END IF
          END IF
          IF(DSMC%ElectronicModel) THEN
            History(iElem,iSpecies,HistTime)%EElec       = -1.0
          ENDIF
        ENDIF
        Sampler(iElem,iSpecies)%Energy(1:3)   = 0
        Sampler(iElem,iSpecies)%Velocity(1:3) = 0
        Sampler(iElem,iSpecies)%PartNum       = 0
        Sampler(iElem,iSpecies)%ERot          = 0
        Sampler(iElem,iSpecies)%EVib          = 0
        Sampler(iElem,iSpecies)%EElec         = 0
      ENDDO
    ENDDO
    iSamplingIters = 0
    IF(HistTime.EQ.nTime) THEN
! Analyze History of averaged properties and check for Steady State Conditions
!  Check all species+properties combinations separately
      Check(1:10) = 0
      DO iSpecies = 1, nSpecies
        ! Check Velocities
        DO iVal = 1, 3
          IF(CheckVelo(iVal)) THEN
            ! Write Velocity Data to Array CheckHistory
            CheckHistory(1:nElems,1:nTime) = History(1:nElems,iSpecies,1:nTime)%Velocity(iVal)
            ! Call Detection Algorithm for Species iSpecies, Property iVal
            CALL SteadyStateDetection_Algorithm(iSpecies,iVal)
            IF(SteadyIdentGlobal(iSpecies,iVal).EQ.1) THEN
              Check(iVal) = Check(iVal)+1
            ENDIF
          ELSE
            Check(iVal) = Check(iVal)+1
          ENDIF
        ENDDO
        ! Check Energies (translational)
        DO iVal = 1, 3
          IF(CheckEnerg(iVal)) THEN
            CheckHistory(1:nElems,1:nTime) = History(1:nElems,iSpecies,1:nTime)%Energy(iVal)
            CALL SteadyStateDetection_Algorithm(iSpecies,iVal+3)
            IF(SteadyIdentGlobal(iSpecies,iVal+3).EQ.1) THEN
              Check(iVal+3) = Check(iVal+3)+1
            ENDIF
          ELSE
            Check(iVal+3) = Check(iVal+3)+1
          ENDIF
        ENDDO
        ! Check Density
        IF(CheckDens) THEN
          CheckHistory(1:nElems,1:nTime) = History(1:nElems,iSpecies,1:nTime)%PartNum
          CALL SteadyStateDetection_Algorithm(iSpecies,7)
          IF(SteadyIdentGlobal(iSpecies,7).EQ.1) THEN
            Check(7) = Check(7)+1
          ENDIF
        ELSE
          Check(7) = Check(7)+1
        ENDIF
        ! Check Energies (internal)
        IF((CollisMode.EQ.2).OR.(CollisMode.EQ.3)) THEN
          IF ((SpecDSMC(iSpecies)%InterID.EQ.2).OR.(SpecDSMC(iSpecies)%InterID.EQ.20)) THEN
            IF(CheckIntEn(1)) THEN
              CheckHistory(1:nElems,1:nTime) = History(1:nElems,iSpecies,1:nTime)%ERot
              CALL SteadyStateDetection_Algorithm(iSpecies,8)
              IF(SteadyIdentGlobal(iSpecies,8).EQ.1) THEN
                Check(8) = Check(8)+1
              ENDIF
            ELSE
              Check(8) = Check(8)+1
            ENDIF
            IF(CheckIntEn(2)) THEN
              CheckHistory(1:nElems,1:nTime) = History(1:nElems,iSpecies,1:nTime)%EVib
              CALL SteadyStateDetection_Algorithm(iSpecies,9)
              IF(SteadyIdentGlobal(iSpecies,9).EQ.1) THEN
                Check(9) = Check(9)+1
              ENDIF
            ELSE
              Check(9) = Check(9)+1
            ENDIF
          ELSE
            Check(8:9) = Check(8:9)+1
          END IF
        ELSE
          Check(8:9) = Check(8:9)+1
        ENDIF
        IF(DSMC%ElectronicModel.AND.CheckIntEn(3)) THEN
          CheckHistory(1:nElems,1:nTime) = History(1:nElems,iSpecies,1:nTime)%EElec
          CALL SteadyStateDetection_Algorithm(iSpecies,10)
          IF(SteadyIdentGlobal(iSpecies,10).EQ.1) THEN
            Check(10) = Check(10)+1
          ENDIF
        ELSE
          Check(10) = Check(10)+1
        ENDIF
      ENDDO
      TotCheck = 0
      DO iVal = 1, 10
        IF(Check(iVal).EQ.nSpecies) THEN
          TotCheck = TotCheck + 1
        ENDIF
      ENDDO
      IF(Totcheck.EQ.10) THEN
        SamplingActive = .TRUE.
      ENDIF
    ENDIF
  ENDIF

END SUBROUTINE SteadyStateDetection_main


SUBROUTINE SteadyStateDetection_Sampling()

  USE MOD_Mesh_Vars,             ONLY : nElems
  USE MOD_Particle_Vars,         ONLY : PEM, usevMPF, PartSpecies, PartMPF, PartState
  USE MOD_DSMC_Vars,             ONLY : Sampler, DSMC, PartStateIntEn, CollisMode, SpecDSMC

!--------------------------------------------------------------------------------------------------------!
! Sampling of macroscopic properties (over nSamplingIters) for reduction of variance and autocorrelation effects
!--------------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                         !
!--------------------------------------------------------------------------------------------------------!
! argument list declaration                                                                              !
! Local variable declaration
  INTEGER           :: iElem, iParts, PartIndex, nParts
!--------------------------------------------------------------------------------------------------------!

  DO iElem = 1, nElems
    nParts = PEM%pNumber(iElem)
    PartIndex = PEM%pStart(iElem)
    DO iParts = 1, nParts
      ! Summation of energy and velocity for each cell and species
      IF(usevMPF) THEN
        Sampler(iElem,PartSpecies(PartIndex))%Energy(1:3)   = Sampler(iElem,PartSpecies(PartIndex))%Energy(1:3) &
                                                            + 0.5 * PartMPF(PartIndex) * PartState(PartIndex,4:6)**2
        Sampler(iElem,PartSpecies(PartIndex))%Velocity(1:3) = Sampler(iElem,PartSpecies(PartIndex))%Velocity(1:3) &
                                                            + PartMPF(PartIndex) * PartState(PartIndex,4:6)
        Sampler(iElem,PartSpecies(PartIndex))%PartNum       = Sampler(iElem,PartSpecies(PartIndex))%PartNum &
                                                            + PartMPF(PartIndex)
        IF((CollisMode.EQ.2).OR.(CollisMode.EQ.3))THEN
          IF((SpecDSMC(PartSpecies(PartIndex))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(PartIndex))%InterID.EQ.20)) THEN
            Sampler(iElem,PartSpecies(PartIndex))%EVib        = Sampler(iElem,PartSpecies(PartIndex))%EVib &
                                                              + PartStateIntEn(PartIndex,1) * PartMPF(PartIndex)
            Sampler(iElem,PartSpecies(PartIndex))%ERot        = Sampler(iElem,PartSpecies(PartIndex))%ERot &
                                                              + PartStateIntEn(PartIndex,2) * PartMPF(PartIndex)
          END IF
        END IF
        IF(DSMC%ElectronicModel) THEN
          Sampler(iElem,PartSpecies(PartIndex))%EElec       = Sampler(iElem,PartSpecies(PartIndex))%EElec &
                                                            + PartStateIntEn(PartIndex,2) * PartMPF(PartIndex)
        ENDIF
      ELSE
        Sampler(iElem,PartSpecies(PartIndex))%Energy(1:3)   = Sampler(iElem,PartSpecies(PartIndex))%Energy(1:3) &
                                                            + 0.5 * PartState(PartIndex,4:6)**2
        Sampler(iElem,PartSpecies(PartIndex))%Velocity(1:3) = Sampler(iElem,PartSpecies(PartIndex))%Velocity(1:3) &
                                                            + PartState(PartIndex,4:6)
        Sampler(iElem,PartSpecies(PartIndex))%PartNum       = Sampler(iElem,PartSpecies(PartIndex))%PartNum &
                                                            + 1
        IF((CollisMode.EQ.2).OR.(CollisMode.EQ.3))THEN
          IF((SpecDSMC(PartSpecies(PartIndex))%InterID.EQ.2).OR.(SpecDSMC(PartSpecies(PartIndex))%InterID.EQ.20)) THEN
            Sampler(iElem,PartSpecies(PartIndex))%EVib        = Sampler(iElem,PartSpecies(PartIndex))%EVib &
                                                              + PartStateIntEn(PartIndex,1)
            Sampler(iElem,PartSpecies(PartIndex))%ERot        = Sampler(iElem,PartSpecies(PartIndex))%ERot &
                                                              + PartStateIntEn(PartIndex,2)
          END IF
        END IF
        IF(DSMC%ElectronicModel) THEN
          Sampler(iElem,PartSpecies(PartIndex))%EElec       = Sampler(iElem,PartSpecies(PartIndex))%EElec &
                                                            + PartStateIntEn(PartIndex,2)
        ENDIF
      ENDIF
      PartIndex = PEM%pNext(PartIndex)
    ENDDO
  ENDDO

END SUBROUTINE SteadyStateDetection_Sampling


SUBROUTINE SteadyStateDetection_Algorithm(iSpec,iVal)

  USE MOD_Mesh_Vars,             ONLY : nElems
  USE MOD_TimeDisc_Vars,         ONLY : iter ! debug
  USE MOD_DSMC_Vars,             ONLY : nTime, CheckHistory!, Sampler, History
  USE MOD_DSMC_Vars,             ONLY : SteadyIdent, SteadyIdentGlobal
  USE MOD_DSMC_Vars,             ONLY : Critical, RValue     ! Von-Neumann-Test
  USE MOD_DSMC_Vars,             ONLY : Epsilon1, Epsilon2, ED_Delta   ! Euclidian Distance Test
  USE MOD_DSMC_Vars,             ONLY : StudCrit, Stud_Indicator    ! Student-t-Test
  USE MOD_DSMC_Vars,             ONLY : PITCrit, ConvCoeff, PIT_Drift   ! Polynomial Interpolation Test
  USE MOD_DSMC_Vars,             ONLY : MK_Trend       ! Mann Kendall (Trend) Test
#ifdef MPI
  USE mpi
#endif

!--------------------------------------------------------------------------------------------------------!
! Steady State Detection based on macroscopic properties using a statistical detection algorithm
!  Identifier 'Method' decides which statistical test to use:
!   1: Von-Neumann-Test
!   2: Euclidian Distance Test
!   3: Student-t-Test
!   4: Polynomial Interpolation Test
!   5: Mann Kendall (Trend) Test
!--------------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                         !
!--------------------------------------------------------------------------------------------------------!
! argument list declaration
  INTEGER, INTENT(IN) :: iSpec  ! Species
  INTEGER, INTENT(IN) :: iVal  ! Property
! Local variable declaration
  INTEGER           :: Method    ! Set Test Method to use
  REAL              :: Alpha    ! Set Level of Significance
  INTEGER, ALLOCATABLE :: NoParts(:)  ! Identifier if there were no particles in cell (number of cells)
  INTEGER           :: NoPartsTest

  REAL              :: Mean      ! Arithmetic Mean of Sample
  REAL              :: Drift      ! Trend estimated via Linear Regression
  REAL, ALLOCATABLE :: MeanDetrended(:)    ! Arithmetic Mean for Data reduced by Trend (number of cells)
  REAL, ALLOCATABLE :: SigmaDetrended(:)  ! Trend Corrected Standard Deviation (sqrt(Variance)) (number of cells)
  REAL,ALLOCATABLE  :: RelSigma(:)    ! Trend Corrected Variation Coefficient (Sigma/Mean) (number of cells)

  ! 1: Von-Neumann-Test
  REAL              :: sVar2, dVar2    ! Sample Variance, Mean Squared Successive Difference
  ! 2: Euclidian Distance Test
  REAL              :: ED       ! Euclidian Distance - Statistic
  REAL              :: VectorLength, VectorMean  ! Euclidian Norm, Arithmetic Norm of the Sample
  REAL              :: ED_Epsilon    ! Critical Value for the ED - Statistic
  ! 3: Student-T-Test
  REAL              :: Indicator    ! Number of Elements of the Sample with Values voted in steady regime
  ! 4: Polynomial Interpolation Test
  REAL              :: DriftP      ! Trend estimated via Savitzky-Golay-Filter
  ! 5: Mann Kendall Test
  REAL              :: Trend      ! Normalized Trend estimated from Mann Kendall Test
  REAL              :: MKCrit      ! Critical Value for the Mann Kendall Test (Normal Distribution)

  INTEGER           :: iCounter, iCounter_global  ! Counter for Cells in steady state (Single, MPI)
  INTEGER           :: nElems_global      ! Overall Number of Cells (MPI)

  INTEGER           :: iElem, iTime, jTime
#ifdef MPI
  INTEGER           :: IERROR
#endif
  INTEGER           :: iProc ! DEBUG
!--------------------------------------------------------------------------------------------------------!

#ifdef MPI
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,iProc,IERROR)
#else
   iProc = 0
#endif
   !WRITE(string,*)iProc


   ! Select Test Method
   Method = 5

   ! Select Test Significance
   !  I.E.: If (nCells_steady/nCells) > (1-Alpha), then accept steady state
   Alpha = 0.02

   IF(iProc.EQ.0) OPEN(UNIT=100,FILE='Neumann.txt',POSITION='APPEND')
   IF(iProc.EQ.0) OPEN(UNIT=102,FILE='ED.txt',POSITION='APPEND')
   IF(iProc.EQ.0) OPEN(UNIT=105,FILE='Stud.txt',POSITION='APPEND')
   IF(iProc.EQ.0) OPEN(UNIT=107,FILE='PIT.txt',POSITION='APPEND')
   IF(iProc.EQ.0) OPEN(UNIT=109,FILE='MK.txt',POSITION='APPEND')

   ALLOCATE(MeanDetrended(nElems))
   ALLOCATE(SigmaDetrended(nElems))
   ALLOCATE(RelSigma(nElems))

   MeanDetrended(1:nElems) = 0
   SigmaDetrended(1:nElems) = 0
   RelSigma(1:nElems) = 0

   ALLOCATE(NoParts(nElems))

   NoParts(1:nElems) = 0

! Check History for Steady State

! Results for each cell are only valid if there were particles in the cell
   DO iElem = 1, nElems
     NoPartsTest = 0
     DO iTime = 1, nTime
       IF(CheckHistory(iElem,iTime).EQ.-1) THEN
         NoPartsTest = 1
       ELSE
         NoPartsTest = 0
         EXIT
       ENDIF
     ENDDO
     IF(NoPartsTest.EQ.1) NoParts(iElem) = 1
   ENDDO

 ! Variance Estimation (Trend Corrected - Trend Estimation via Linear Regression) - for Methods 2,3,4
   IF((Method.EQ.2).OR.(Method.EQ.3).OR.(Method.EQ.4)) THEN
      DO iElem = 1, nElems

        ! 1) Arithmetic Mean of Sample:
        Mean  = 0
        DO iTime = 1, nTime
          Mean = Mean + CheckHistory(iElem,iTime)
        ENDDO
        Mean = Mean / nTime

        ! 2) Trend Estimation via Linear Regression:
        Drift = 0
        DO iTime = 1, nTime
          Drift = Drift + (iTime - 0.5*(nTime+1)) * (CheckHistory(iElem,iTime) - Mean)
        ENDDO
        Drift = Drift * 12.0/((nTime-1)*nTime*(nTime+1))

        ! 3) Arithmetic Mean for Data reduced by Trend
        DO iTime = 1, nTime
          MeanDetrended(iElem) = MeanDetrended(iElem) + (CheckHistory(iElem,iTime) - Drift*(iTime-1))
        ENDDO
        MeanDetrended(iElem) = MeanDetrended(iElem) / nTime

        ! 4) Trend Corrected Standard Deviation and Variation Coefficient
        DO iTime = 1, nTime
          SigmaDetrended(iElem) = SigmaDetrended(iElem) + (CheckHistory(iElem,iTime) - MeanDetrended(iElem) - Drift*(iTime-1))**2
        ENDDO
        SigmaDetrended(iElem) = sqrt(SigmaDetrended(iElem)/(nTime-2))
        RelSigma(iElem) = abs(SigmaDetrended(iElem)/Mean)

      ENDDO
   ENDIF
 ! End Variance Estimation

   IF(Method.EQ.1) THEN
 ! Steady State Detection with the Von Neumann - Method
 !  Trend detection by analyzing the ratio of the sample variance 'sVar2' and the mean squared successive difference 'dVar2'
 !  Von-Neumann-Ratio: R = 0.5 * s^2 / d^2
 ! Check each Cell separately
     DO iElem = 1, nElems
       Mean = 0
       sVar2 = 0
       dVar2 = 0
       RValue(iElem) = -1.0
       DO iTime = 1, nTime
         Mean = Mean + CheckHistory(iElem,iTime)
       ENDDO
       Mean = Mean/ nTime
       DO iTime = 1, nTime-1
         sVar2 = sVar2 + ( CheckHistory(iElem,iTime) - Mean )**2
         dVar2 = dVar2 + ( CheckHistory(iElem,iTime+1) - CheckHistory(iElem,iTime) )**2
       ENDDO
       sVar2 = sVar2 + ( CheckHistory(iElem,nTime) - Mean )**2
       sVar2 = sVar2 / (nTime - 1)
       dVar2 = dVar2 / (nTime - 1)
       RValue(iElem) = 0.5 * dVar2 / sVar2
       IF(iProc.EQ.0.AND.iElem.EQ.1) WRITE(100,*)iter+1,RValue(iElem)
       IF((RValue(iElem).GT.Critical(1)).AND.(RValue(iElem).LT.Critical(2))) THEN
         SteadyIdent(iElem,iSpec,iVal) = 1
       ELSE
         SteadyIdent(iElem,iSpec,iVal) = 0
       ENDIF
     ENDDO
   ELSEIF(Method.EQ.2) THEN
 ! Steady State Detection with Euclidian Distance Method
 !  Trend Detection based on the Quotient of the Arithmetic Mean and the Euclidian Norm of the Sample
 ! Check each Cell separately
     DO iElem = 1, nElems
       VectorLength = 0
       VectorMean = 0
       ED = 0
       ED_Epsilon = Epsilon1 * (RelSigma(iElem)**Epsilon2)
       DO iTime = 1, nTime
         VectorLength = VectorLength + CheckHistory(iElem,iTime)**2
         VectorMean   = VectorMean + CheckHistory(iElem,iTime)
       ENDDO
       VectorLength = sqrt(VectorLength)
       VectorMean = VectorMean / nTime
       ED = abs( VectorMean/VectorLength * sqrt(REAL(nTime)) )
       IF(iProc.EQ.0.AND.iElem.EQ.1) WRITE(102,*)iter+1,ED_Epsilon,ED
       ED_Delta(iElem) = (1 - ED) / ED_Epsilon
       IF((ED.GT.(1.0-ED_Epsilon)).AND.(ED.LT.(1.0+ED_Epsilon))) THEN
         SteadyIdent(iElem,iSpec,iVal) = 1
       ELSE
         SteadyIdent(iElem,iSpec,iVal) = 0
       ENDIF
     ENDDO
   ELSEIF(Method.EQ.3) THEN
 ! Steady State Detection with Student-t Test
 !  Trend Detection with Student-t-Test on the difference between actual value and mean with trend corrected data
 ! Check each Cell separately
     DO iElem = 1, nElems
       Indicator = 0
       DO iTime = 1, nTime
         IF( abs(CheckHistory(iElem,iTime)-MeanDetrended(iElem)).LE.(StudCrit*SigmaDetrended(iElem)) ) THEN
           Indicator = Indicator + 1
         ENDIF
       ENDDO
       Indicator = Indicator / nTime
       IF(iProc.EQ.0.AND.iElem.EQ.1) WRITE(105,*)iter+1,Indicator
       Stud_Indicator(iElem) = Indicator
       IF(Indicator.GE.0.99) THEN
         SteadyIdent(iElem,iSpec,iVal) = 1
       ELSE
         SteadyIdent(iElem,iSpec,iVal) = 0
       ENDIF
     ENDDO
   ELSEIF(Method.EQ.4) THEN
! Steady State Detection with Polynomial Interpolation Test (based on Savitzky & Golay Filter)
!  Trend computed as Slope of a Cubic Least Square Interpolation Polynom through nTime Points
! Check each Cell separately
     DO iElem = 1, nElems
       DriftP = 0
       DO iTime = 1, nTime
         DriftP = DriftP + CheckHistory(iElem,iTime) * ConvCoeff(iTime)
       ENDDO
       IF(abs(DriftP).LT.(SigmaDetrended(iElem)*PITCrit)) THEN
         SteadyIdent(iElem,iSpec,iVal) = 1
       ELSE
         SteadyIdent(iElem,iSpec,iVal) = 0
       ENDIF
       IF(iProc.EQ.0.AND.iElem.EQ.1) WRITE(107,*)iter+1,abs(DriftP),SigmaDetrended(iElem)*PITCrit
       PIT_Drift(iElem) = abs(DriftP) / ( SigmaDetrended(iElem) * PITCrit )
     ENDDO
   ELSEIF(Method.EQ.5) THEN
 ! Steady State Detection with the Mann Kendall Test
 !  Rank-based Test on Significance of Trend
 ! Check each Cell separately
     DO iElem = 1, nElems
       Trend = 0
       DO iTime = 1, nTime-1
         DO jTime = iTime+1, nTime
           IF((CheckHistory(iElem,jTime)-CheckHistory(iElem,iTime)).GT.0) THEN
             Trend = Trend + 1
           ELSEIF((CheckHistory(iElem,jTime)-CheckHistory(iElem,iTime)).LT.0) THEN
             Trend = Trend - 1
           ENDIF
         ENDDO
       ENDDO
       IF(Trend.GT.0) THEN
         Trend = (Trend - 1)/sqrt(REAL(nTime*(nTime-1)*(2*nTime+5))/18)
       ELSEIF(Trend.LT.0) THEN
         Trend = (Trend + 1)/sqrt(REAL(nTime*(nTime-1)*(2*nTime+5))/18)
       ENDIF
       !MKCrit = 1.96    ! 5% Significance
       MKCrit = 2.3263    ! 2% Significance
       !MKCrit = 2.5758    ! 1% Significance
       MK_Trend(iElem) = Trend/MKCrit
       IF(iProc.EQ.0.AND.iElem.EQ.1) WRITE(109,*)iter+1,MK_Trend(iElem)
       IF(abs(Trend).LT.MKCrit) THEN
         SteadyIdent(iElem,iSpec,iVal) = 1
       ELSE
         SteadyIdent(iElem,iSpec,iVal) = 0
       ENDIF
     ENDDO
   ENDIF

 ! Check if all Cells are in Steady State
   iCounter = 0
   DO iElem = 1, nElems
     IF((SteadyIdent(iElem,iSpec,iVal).EQ.1).AND.(NoParts(iElem).EQ.0)) iCounter = iCounter + 1
   ENDDO
   !
#ifdef MPI
   CALL MPI_ALLREDUCE(iCounter, iCounter_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERROR)
   CALL MPI_ALLREDUCE(nElems, nElems_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, IERROR)
#else
   iCounter_global = iCounter
   nElems_global = nElems
#endif
   !
   IF(iCounter_global.GE.(nElems_global*(1.0-Alpha))) SteadyIdentGlobal(iSpec,iVal) = 1

   DEALLOCATE(MeanDetrended)
   DEALLOCATE(SigmaDetrended)
   DEALLOCATE(RelSigma)

   IF(iProc.EQ.0) CLOSE(100)
   IF(iProc.EQ.0) CLOSE(102)
   IF(iProc.EQ.0) CLOSE(105)
   IF(iProc.EQ.0) CLOSE(107)
   IF(iProc.EQ.0) CLOSE(109)

END SUBROUTINE SteadyStateDetection_Algorithm


SUBROUTINE  EntropyCalculation()

  USE MOD_Mesh_Vars,             ONLY : nElems
  USE MOD_Particle_Vars,         ONLY : PEM, usevMPF, PartState!, Species, PartSpecies, PartMPF
  USE MOD_DSMC_Vars,             ONLY : HValue

!--------------------------------------------------------------------------------------------------!
! Cellwise Calculation of Entropy Parameter 'H' using in Boltzmann's H - Theorem:
!  H = - Integral(f(v)*ln(f(v)),dv,R^3)
!  1. Estimate the Probability Density Function f(v) for the Particle Velocities
!  2. Solve the Integral for H numerically -> HValue
!   a. with the Histogram     (Mode = 1)
!   b. with Kernel Density Estimation   (Mode = 2)
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER               :: Mode  ! Density Estimation Method, 1 = Histogram,  2 = Kernel Density Estimation
  REAL                  :: minvel(3), maxvel(3), avvel(3)  ! Minimum / Maximum / Average Particle Velocity (x,y,z)
  REAL                  :: sigmavel(3)        ! Standard Deviation of Particle Velocity (x,y,z)
  INTEGER               :: iElem, iPart, nParts, PartIndex
 ! Histogram
  INTEGER               :: nBins(3)        ! Number of Intervals in each direction (x,y,z)
  REAL                  :: BinWidth(3)        ! Interval/Bin Size (x,y,z)
  INTEGER, ALLOCATABLE  :: BinCounter(:,:,:)      ! Counter for Number of Particles in each Bin (index x,y,z)
  INTEGER               :: BinIdent(3), iBin, jBin, kBin
 ! Kernel Density Estimation
  REAL                  :: Bandwidth2(3)      ! Bandwidth for Kernel (x,y,z)
  INTEGER               :: nSteps        ! Number of Steps for numerical integration (same for each direction)
  REAL                  :: deltav(3)        ! Step Size for numerical integration (x,y,z)
  REAL                  :: vpos(3)        ! local position in velocity space (x,y,z)
  REAL                  :: vdens        ! local estimate of the density function
  INTEGER               :: iStep, jStep, kStep
!--------------------------------------------------------------------------------------------------!

! Set the Density Estimation Method
  Mode = 1

  ! Loop over all Cells
  DO iElem = 1, nElems

    nParts = PEM%pNumber(iElem)
    IF(nParts.LT.2) THEN
      ! Not enough particles in Cell, Set H to Zero
      HValue(iElem) = 0.0
      CYCLE
    ENDIF

    ! Determine Minimum, Maximum and Average Particle Velocities
    PartIndex = PEM%pStart(iElem)
    maxvel(1:3) = PartState(PartIndex,4:6)
    minvel(1:3) = PartState(PartIndex,4:6)
    avvel(1:3)  = PartState(PartIndex,4:6)
    iPart = 1
    PartIndex = PEM%pNext(PartIndex)
    DO iPart = 2, nParts
      IF(PartState(PartIndex,4).GT.maxvel(1)) maxvel(1) = PartState(PartIndex,4)
      IF(PartState(PartIndex,4).LT.minvel(1)) minvel(1) = PartState(PartIndex,4)
      IF(PartState(PartIndex,5).GT.maxvel(2)) maxvel(2) = PartState(PartIndex,5)
      IF(PartState(PartIndex,5).LT.minvel(2)) minvel(2) = PartState(PartIndex,5)
      IF(PartState(PartIndex,6).GT.maxvel(3)) maxvel(3) = PartState(PartIndex,6)
      IF(PartState(PartIndex,6).LT.minvel(3)) minvel(3) = PartState(PartIndex,6)
      avvel(1:3)  = avvel(1:3) + PartState(PartIndex,4:6)
      PartIndex = PEM%pNext(PartIndex)
    ENDDO
    avvel(1:3)  = avvel(1:3) / nParts

    ! Calculate the Standard Deviation of the Particle Velocities
    ! (for Determination of Optimum Bin-Size resp. Bandwidth)
    PartIndex = PEM%pStart(iElem)
    sigmavel(1:3) = (PartState(PartIndex,4:6) - avvel(1:3))**2
    iPart = 1
    PartIndex = PEM%pNext(PartIndex)
    DO iPart = 2, nParts
      sigmavel(1:3) = sigmavel(1:3) + (PartState(PartIndex,4:6) - avvel(1:3))**2
      PartIndex = PEM%pNext(PartIndex)
    ENDDO
    sigmavel(1:3) = sqrt(sigmavel(1:3)/(nParts-1))

    IF(Mode.EQ.1) THEN
      ! Histogram: Adapt the Velocity Intervals to ensure all paticles are captured
      minvel(1:3) = minvel(1:3) - abs(0.01*minvel(1:3))
      maxvel(1:3) = maxvel(1:3) + abs(0.01*maxvel(1:3))
    ELSE
      ! KDE: Widen the Velocity Intervals to ensure correct integration of the density function
      minvel(1:3) = minvel(1:3) - 2.6 * sigmavel(1:3)
      maxvel(1:3) = maxvel(1:3) + 2.6 * sigmavel(1:3)
    ENDIF

    IF(Mode.EQ.1) THEN
     ! Use the Histogram Method
     !  The Velocity Space (from mimimum to maximum particle velocity in all directions only, i.e. a Cuboid)
     !  is in each direction divided into 'nBins' intervals with constant size (BinWidth).
     !  Thus, we have nBins(1)*nBins(2)*nBins(3) Bins (Cells).
     !  For each Bin with index i,j,k the Number of Particles with velocities lying inside this Bin
     !  is determined (BinCounter(i,j,k)),
     !  as an estimate for the local mean value of the density function; it is
     !  f(i,j,k) = BinCounter(i,j,k)/(nParts*BinWidth(1)*BinWidth(2)*BinWidth(3)).
     !  Here, the Number of Particles is used directly to solve the integral for H with the (3-dimensional) rectangle method:
     !  H = - Summation(f(i,j,k)*ln(f(i,j,k)*DV)

      ! Set Optimum Bin Size (Scott)
      BinWidth(1:3) = 3.5 * sigmavel(1:3) / REAL(nParts)**0.2
      nBins(1:3) = CEILING((maxvel(1:3) - minvel(1:3)) / BinWidth(1:3))
      BinWidth(1:3) = (maxvel(1:3) - minvel(1:3)) / nBins(1:3)

      ALLOCATE(BinCounter(nBins(1),nBins(2),nBins(3)))

      ! Distribute the Particles to the Bins
      BinCounter(1:nBins(1),1:nBins(2),1:nBins(3)) = 0
      PartIndex = PEM%pStart(iElem)
      DO iPart = 1, nParts
        IF(usevMPF) THEN
          !Species(PartSpecies(PartIndex))%MassIC * PartMPF(PartIndex) * PartState(PartIndex,4)
          DO iBin = 1, nBins(1)
            IF((PartState(PartIndex,4).GE.(minvel(1)+BinWidth(1)*(iBin-1))) &
            .AND.(PartState(PartIndex,4).LT.(minvel(1)+BinWidth(1)*iBin))) THEN
              BinIdent(1) = iBin
              EXIT
            ENDIF
          ENDDO
          DO jBin = 1, nBins(2)
            IF((PartState(PartIndex,5).GE.(minvel(2)+BinWidth(2)*(jBin-1))) &
            .AND.(PartState(PartIndex,5).LT.(minvel(2)+BinWidth(2)*jBin))) THEN
              BinIdent(2) = jBin
              EXIT
            ENDIF
          ENDDO
          DO kBin = 1, nBins(3)
            IF((PartState(PartIndex,6).GE.(minvel(3)+BinWidth(3)*(kBin-1))) &
            .AND.(PartState(PartIndex,6).LT.(minvel(3)+BinWidth(3)*kBin))) THEN
              BinIdent(3) = kBin
              EXIT
            ENDIF
          ENDDO
          BinCounter(BinIdent(1),BinIdent(2),BinIdent(3)) = BinCounter(BinIdent(1),BinIdent(2),BinIdent(3)) + 1
        ELSE
          DO iBin = 1, nBins(1)
            IF((PartState(PartIndex,4).GE.(minvel(1)+BinWidth(1)*(iBin-1))) &
            .AND.(PartState(PartIndex,4).LT.(minvel(1)+BinWidth(1)*iBin))) THEN
              BinIdent(1) = iBin
              EXIT
            ENDIF
          ENDDO
          DO jBin = 1, nBins(2)
            IF((PartState(PartIndex,5).GE.(minvel(2)+BinWidth(2)*(jBin-1))) &
            .AND.(PartState(PartIndex,5).LT.(minvel(2)+BinWidth(2)*jBin))) THEN
              BinIdent(2) = jBin
              EXIT
            ENDIF
          ENDDO
          DO kBin = 1, nBins(3)
            IF((PartState(PartIndex,6).GE.(minvel(3)+BinWidth(3)*(kBin-1))) &
            .AND.(PartState(PartIndex,6).LT.(minvel(3)+BinWidth(3)*kBin))) THEN
              BinIdent(3) = kBin
              EXIT
            ENDIF
          ENDDO
          BinCounter(BinIdent(1),BinIdent(2),BinIdent(3)) = BinCounter(BinIdent(1),BinIdent(2),BinIdent(3)) + 1
        ENDIF
        PartIndex = PEM%pNext(PartIndex)
      ENDDO

      ! Calculate the H - Integral with the (3-dimensional) rectangle method
      HValue(iElem) = 0.0
      DO iBin = 1, nBins(1)
        DO jBin = 1, nBins(2)
          DO kBin = 1, nBins(3)
            IF(BinCounter(iBin,jBin,kBin).GT.0) HValue(iElem) = HValue(iElem) - REAL(BinCounter(iBin,jBin,kBin))/nParts &
               * log(REAL(BinCounter(iBin,jBin,kBin))/(nParts*BinWidth(1)*BinWidth(2)*BinWidth(3)))
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(BinCounter)

    ELSE
     ! Use Kernel Density Estimation

      ! Set Optimum Bandwidth for KDE (Silverman)
      Bandwidth2(1:3) = 0.9382 * (sigmavel(1:3)**2) / (REAL(nParts)**0.2857)

      ! Set Stepsize and Boundaries for numerical integration
      nSteps = 12
      deltav(1:3) = (maxvel(1:3) - minvel(1:3))/nSteps

      ! Solve integral numerically using (3-dimensional) rectangle method
      HValue(iElem) = 0.0
      vpos(1) = minvel(1) + 0.5 * deltav(1)
      DO iStep = 1, nSteps
        vpos(2) = minvel(2) + 0.5 * deltav(2)
        DO jStep = 1, nSteps
          vpos(3) = minvel(3) + 0.5 * deltav(3)
          DO kStep = 1, nSteps
            ! Caluclate the Value of the Density Function (vdens) at the actual Position (vpos) in Velocity Space
            CALL KernelDensityEstimation(iElem,vpos,Bandwidth2,vdens)
            IF(vdens.NE.0) HValue(iElem) = HValue(iElem) - vdens * log(vdens) * deltav(1) * deltav(2) * deltav(3)
            vpos(3) = vpos(3) + deltav(3)
          ENDDO
          vpos(2) = vpos(2) + deltav(2)
        ENDDO
        vpos(1) = vpos(1) + deltav(1)
      ENDDO

    ENDIF

  ENDDO

END SUBROUTINE EntropyCalculation

SUBROUTINE KernelDensityEstimation(Elem,Pos,Bandwidth2,Dens)

  USE MOD_Particle_Vars,         ONLY : PEM, PartState

!--------------------------------------------------------------------------------------------------!
! Kernel Density Estimation (3D Gaussian Kernel)
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration
  INTEGER, INTENT(IN)   :: Elem
  REAL, INTENT(IN)      :: Pos(3)    ! Position (x,y,z)
  REAL, INTENT(IN)      :: Bandwidth2(3)  ! Kernel Bandwidth (x,y,z)
  REAL, INTENT(OUT)     :: Dens      ! Estimated local Density
! Local variable declaration                                                                       !
  INTEGER               :: iPart, nParts, PartIndex
  REAL                  :: PI
!--------------------------------------------------------------------------------------------------!

  PI = 4*ATAN(1.0)
  Dens = 0.0
  nParts = PEM%pNumber(Elem)
  PartIndex = PEM%pStart(Elem)

  DO iPart = 1, nParts
    Dens = Dens + exp(-0.5*(((Pos(1)-PartState(PartIndex,4))**2)/Bandwidth2(1)+((Pos(2)-PartState(PartIndex,5))**2) &
         / Bandwidth2(2)+((Pos(3)-PartState(PartIndex,6))**2)/Bandwidth2(3)))
    PartIndex = PEM%pNext(PartIndex)
  ENDDO
  Dens = Dens / (nParts*((2*PI)**1.5)*sqrt(Bandwidth2(1)*Bandwidth2(2)*Bandwidth2(3)))

END SUBROUTINE KernelDensityEstimation


END MODULE MOD_DSMC_SteadyState
