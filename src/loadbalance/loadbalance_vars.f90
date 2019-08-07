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
MODULE MOD_LoadBalance_Vars
!===================================================================================================================================
! Variables needed for the evaluation of the record points
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                             :: DoLoadBalance                              ! DoLoadBalance
LOGICAL                             :: PerformLoadBalance=.FALSE.                 ! Flag if loadbalance is performed in current iter
INTEGER                             :: LoadBalanceSample                          ! Number of samples for loadbalance
LOGICAL                             :: PerformLBSample                            ! Flag for enabling time measurement in current
                                                                                  ! timestep (automatically set depending on LB
                                                                                  ! sampling method)
LOGICAL                             :: PerformPartWeightLB                        ! Flag for performing LB with partMPIWeight
                                                                                  ! instead of summed Elemtimes
                                                                                  ! -> nParts*PartWeight written into elemtime array
LOGICAL                             :: InitLoadBalanceIsDone                      ! switch for checking

! time measurement
REAL,ALLOCATABLE                    :: tCurrent(:)                                ! time measurement over one step
                                                                                  ! measured elem-independent and later weighted
                                                                                  ! for indeces look into piclas.h
REAL,ALLOCATABLE                    :: tCurrent_LB_DG(:)                                ! time measurement over one step
! counter
INTEGER                             :: nLoadBalance                               ! number of load balances
INTEGER                             :: nLoadBalanceSteps                          ! number of performed  load balances steps
REAL,ALLOCATABLE                    :: LoadDistri(:)                              ! Weighted load distribution of all procs
INTEGER,ALLOCATABLE                 :: PartDistri(:)                              ! Part distribution of all procs
REAL                                :: MaxWeight                                  ! Maximum Weight of proc on domain
REAL                                :: MinWeight                                  ! Minimum Weight of proc on domain
REAL                                :: CurrentImbalance
REAL                                :: NewImbalance                               ! Imbalance after rebalance step

TYPE tData
  INTEGER, ALLOCATABLE :: offsetElemMPI(:)
  INTEGER              :: numOfCalls
  TYPE(tData), POINTER :: nextData => null()
END TYPE tData
TYPE(tData), POINTER :: firstData => null() !linked-list of old offsetElemMPI for WeightDistributionMethod 5 and 6

!-----------------------------------------------------------------------------------------------------------------------------------
! particle load balancing
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER(KIND=8)                     :: nSkipAnalyze                               ! Skip Analyze-Dt
REAL                                :: ParticleMPIWeight
REAL                                :: DeviationThreshold                         ! threshold for load-balancing
LOGICAL                             :: writePartitionInfo                         ! write partitioninfo file
REAL                                :: WeightSum                                  ! global sum of all weights
REAL                                :: targetWeight                               ! optimal weight for each proc

!-----------------------------------------------------------------------------------------------------------------------------------
! Element Local measurement
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                    :: ElemTime(:)
INTEGER,ALLOCATABLE                 :: ElemHDGSides(:) ! number of master sides for the HDG solver for each element
INTEGER                             :: TotalHDGSides   ! total number of master sides for the HDG solver over all local elements
REAL,ALLOCATABLE                    :: ElemGlobalTime(:)
INTEGER(KIND=8),ALLOCATABLE         :: nPartsPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nDeposPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nTracksPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nSurfacefluxPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nPartsPerBCElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nSurfacePartsPerElem(:)


END MODULE MOD_LoadBalance_Vars
