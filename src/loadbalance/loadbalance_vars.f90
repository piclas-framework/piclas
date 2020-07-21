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
LOGICAL                             :: DoLoadBalance               ! DoLoadBalance
LOGICAL                             :: PerformLoadBalance=.FALSE.  ! Flag if load balance is performed in current time step iteration
INTEGER                             :: LoadBalanceSample           ! Number of samples for load balance
LOGICAL                             :: PerformLBSample             ! Flag for enabling time measurement in current
                                                                   ! Time step (automatically set depending on LB
                                                                   ! sampling method)
LOGICAL                             :: PerformPartWeightLB         ! Flag for performing LB with partMPIWeight
                                                                   ! instead of summed ElemTimes
                                                                   ! -> nParts*PartWeight written into elemtime array
LOGICAL                             :: InitLoadBalanceIsDone       ! Switch for checking

! time measurement
REAL,ALLOCATABLE                    :: tCurrent(:)                 ! Time measurement over one step
                                                                   ! measured elem-independent and later weighted
                                                                   ! for indices look into piclas.h

REAL,ALLOCATABLE                    :: tCurrent_LB_DG(:)           ! Time measurement over one step
! counter
INTEGER                             :: nLoadBalance                ! Number of load balances calculations (calls of ComputeElemLoad)
INTEGER                             :: nLoadBalanceSteps           ! Number of performed load balances steps
INTEGER                             :: LoadBalanceMaxSteps         ! Number of maximum allowed performed load balances steps
REAL,ALLOCATABLE                    :: LoadDistri(:)               ! Weighted load distribution of all procs
INTEGER,ALLOCATABLE                 :: PartDistri(:)               ! Part distribution of all procs
REAL                                :: MaxWeight                   ! Maximum Weight of proc on domain
REAL                                :: MinWeight                   ! Minimum Weight of proc on domain
REAL                                :: CurrentImbalance
REAL                                :: NewImbalance                ! Imbalance after rebalance step

TYPE tData
  INTEGER, ALLOCATABLE :: offsetElemMPI(:)
  INTEGER              :: numOfCalls
  TYPE(tData), POINTER :: nextData => null()
END TYPE tData
TYPE(tData), POINTER :: firstData => null() !linked-list of old offsetElemMPI for WeightDistributionMethod 5 and 6

!-----------------------------------------------------------------------------------------------------------------------------------
! particle load balancing
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER(KIND=8)                     :: nSkipAnalyze                ! Skip Analyze-Dt
REAL                                :: ParticleMPIWeight
REAL                                :: DeviationThreshold          ! threshold for load-balancing
LOGICAL                             :: writePartitionInfo          ! write partitioninfo file
REAL                                :: WeightSum                   ! global sum of all weights
REAL                                :: targetWeight                ! optimal weight for each proc

!-----------------------------------------------------------------------------------------------------------------------------------
! Element Local measurement
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE                    :: ElemTime(:)
LOGICAL                             :: NullifyElemTime
REAL,ALLOCATABLE                    :: ElemTime_tmp(:)  ! Additional container for restarting and keeping the old ElemTime values in
                                                        ! the state.h5 file
REAL                                :: ElemTimePartTot  ! Total time spent for particle routines (all procs)
REAL                                :: ElemTimeFieldTot ! Total time spent for field routines (all procs)
REAL                                :: ElemTimePart     ! Time spent for particle routines
REAL                                :: ElemTimeField    ! Time spent for field routines
REAL,ALLOCATABLE                    :: ElemHDGSides(:)  ! number of master sides for the HDG solver for each element
REAL                                :: TotalHDGSides    ! total number of master sides for the HDG solver over all local elements
REAL,ALLOCATABLE                    :: ElemGlobalTime(:)
INTEGER(KIND=8),ALLOCATABLE         :: nPartsPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nDeposPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nTracksPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nSurfacefluxPerElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nPartsPerBCElem(:)
INTEGER(KIND=8),ALLOCATABLE         :: nSurfacePartsPerElem(:)


END MODULE MOD_LoadBalance_Vars
