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
!===================================================================================================================================
!> Contains global variables used by the Timedisc modules.
!===================================================================================================================================
MODULE MOD_TimeDisc_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL             :: time_start                        !< start time of the timedisc
REAL             :: ManualTimeStep                    !> Manual TimeStep
LOGICAL          :: useManualTimeStep                 !> Logical Flag for manual timestep. For consistency
                                                      !> with IAG programming style
REAL             :: tWallRemaining                    !> Remaining wall time in hours
REAL             :: tStart                            !> simulation time at the beginning of the simulation
REAL             :: tEnd                              !> simulation end time
REAL             :: tAnalyze                          !> time for next analyze
REAL             :: time                              !> Simulation Time
REAL             :: dt                                !> simulation time step
REAL             :: dtWeight                          !> part of original dt that is currently used as dt (output, rk, ...)
REAL             :: CFLScale                          !> cfl scale
REAL             :: CFLtoOne                          !> scaling factor to scale CFL to one
REAL             :: sdtCFLOne                         !> inverse of dt of CFLOne
REAL             :: RKdtFrac,RKdtFracTotal
INTEGER          :: iStage
INTEGER(KIND=8)  :: iter                              !> iteration since first init
INTEGER(KIND=8)  :: IterDisplayStep                   !> number of displayed iteration steps written during simulation
INTEGER(KIND=8)  :: IterDisplayStepUser               !> number of displayed iteration steps that are defined by user
LOGICAL          :: DoDisplayIter                     !> flag if iterations are displayed (TRUE if IterDisplayStep>0)
LOGICAl          :: TimediscInitIsDone = .FALSE.
REAL             :: TimeDG, TimeParticle
#if defined(PARTICLES) && USE_HDG
REAL             :: dt_Min(4) !> dt_Min(DT_MIN)       = dt_Min(1) = original dt_Min
                              !> dt_Min(DT_ANALYZE)   = dt_Min(2) = tAnalyzeDiff
                              !> dt_Min(DT_END)       = dt_Min(3) = tEndDiff
                              !> dt_Min(DT_BR_SWITCH) = dt_Min(4) = tBRDiff (time to BR<->kin switch)
#else
REAL             :: dt_Min(3) !> dt_Min(DT_MIN)       = dt_Min(1) = original dt_Min
                              !> dt_Min(DT_ANALYZE)   = dt_Min(2) = tAnalyzeDiff
                              !> dt_Min(DT_END)       = dt_Min(3) = tEndDiff
#endif /*defined(PARTICLES) && USE_HDG*/
!REAL             :: tEndDiff     !> difference between simulation time and simulation end time -> dt_Min(DT_END)
!REAL             :: tAnalyzeDiff !> difference between simulation time and next analyze time -> dt_Min(DT_ANALYZE)
#if (PP_TimeDiscMethod==509)
REAL             :: dt_old
#endif /*(PP_TimeDiscMethod==509)*/
#if (PP_TimeDiscMethod==100)
INTEGER,PARAMETER:: nRKStages=1
#endif
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)|| (PP_TimeDiscMethod==6)
! DG solution vol
TYPE N_Ut_Vol
  REAL,ALLOCATABLE :: Ut_temp( :,:,:,:) !> Temporal variable for Ut
  REAL,ALLOCATABLE :: U2t_temp(:,:,:,:) !> Temporal variable for U2t (PML)
END TYPE N_Ut_Vol

! DG solution (JU or U) vectors)
TYPE(N_Ut_Vol),ALLOCATABLE :: Ut_N(:)       !< Solution variable for each equation, node and element,

#endif

!-----------------------------------------------------------------------------------------------------------------------------------
! TIME INTEGRATION: RUNGE_KUTTA COEFFICIENTS AND STABILITY NUMBERS
!-----------------------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==501)
! Low-storage Runge-Kutta 3, 3 stages, Kopriva,Algorithm 42
CHARACTER(LEN=255),PARAMETER :: TimeDiscName = 'STANDARD RK3-3'
INTEGER,PARAMETER  :: nRKStages=3
REAL,PARAMETER  :: RK3_a2= 5./9.
REAL,PARAMETER  :: RK3_a3= 153./128.
REAL,PARAMETER  :: RK_a(2:3) = (/RK3_a2,RK3_a3/)

REAL,PARAMETER  :: RK3_b1= 1./3.
REAL,PARAMETER  :: RK3_b2= 15./16.
REAL,PARAMETER  :: RK3_b3= 8./15.
REAL,PARAMETER  :: RK_b(1:3) = (/RK3_b1,RK3_b2,RK3_b3/)

REAL,PARAMETER  :: RK3_c2= 1./3.
REAL,PARAMETER  :: RK3_c3= 0.75
REAL,PARAMETER  :: RK_c(2:3) = (/RK3_c2,RK3_c3/)

! define scaling factor for CFL for degree N=1..15
REAL,PARAMETER  :: CFLScaleAlpha(1:15) = &
#if (PP_NodeType==1)
(/ 1.2285, 1.0485, 0.9101, 0.8066, 0.7268, 0.6626, 0.6109, 0.5670, 0.5299, 0.4973, 0.4703, 0.4455, 0.4230, 0.4039, 0.3859 /)
#elif (PP_NodeType==2)
(/ 3.1871, 2.2444, 1.7797, 1.5075, 1.3230, 1.1857, 1.0800, 0.9945, 0.9247, 0.8651, 0.8134, 0.7695, 0.7301, 0.6952, 0.6649 /)
#endif /*PP_NodeType*/
#endif

#if ((PP_TimeDiscMethod==4))
INTEGER,PARAMETER  :: nRKStages=1
#endif

#if ((PP_TimeDiscMethod==2)||(PP_TimeDiscMethod==502))
! Runge-Kutta 4 - Carpenter 1994 NASA Report
INTEGER,PARAMETER  :: nRKStages=5
REAL,PARAMETER  :: RK4_a2=  567301805773.0/  1357537059087.0
REAL,PARAMETER  :: RK4_a3= 2404267990393.0/  2016746695238.0
REAL,PARAMETER  :: RK4_a4= 3550918686646.0/  2091501179385.0
REAL,PARAMETER  :: RK4_a5= 1275806237668.0/   842570457699.0
REAL,PARAMETER  :: RK_a(2:5) = (/RK4_a2,RK4_a3,RK4_a4,RK4_a5/)

REAL,PARAMETER  :: RK4_b1= 1432997174477.0/  9575080441755.0
REAL,PARAMETER  :: RK4_b2= 5161836677717.0/ 13612068292357.0
REAL,PARAMETER  :: RK4_b3= 1720146321549.0/  2090206949498.0
REAL,PARAMETER  :: RK4_b4= 3134564353537.0/  4481467310338.0
REAL,PARAMETER  :: RK4_b5= 2277821191437.0/ 14882151754819.0
REAL,PARAMETER  :: RK_b(1:5) = (/RK4_b1,RK4_b2,RK4_b3,RK4_b4,RK4_b5 /)

REAL,PARAMETER  :: RK4_c2= 1432997174477.0/  9575080441755.0
REAL,PARAMETER  :: RK4_c3= 2526269341429.0/  6820363962896.0
REAL,PARAMETER  :: RK4_c4= 2006345519317.0/  3224310063776.0
REAL,PARAMETER  :: RK4_c5= 2802321613138.0/  2924317926251.0
REAL,PARAMETER  :: RK_c(2:5) = (/RK4_c2,RK4_c3,RK4_c4,RK4_c5/)

! for degree N=1..15
REAL,PARAMETER  :: CFLScaleAlpha(1:15) = &
#if (PP_NodeType==1)
(/ 2.0351, 1.7595, 1.5401, 1.3702, 1.2375, 1.1318, 1.0440, 0.9709, 0.9079, 0.8539, 0.8066, 0.7650, 0.7290, 0.6952, 0.6660 /)
#elif (PP_NodeType==2)
(/ 4.7497, 3.4144, 2.8451, 2.4739, 2.2027, 1.9912, 1.8225, 1.6830, 1.5682, 1.4692, 1.3849, 1.3106, 1.2454, 1.1880, 1.1362 /)
#endif /*PP_NodeType*/
#endif

#if (PP_TimeDiscMethod==6)||(PP_TimeDiscMethod==506)
! Low storage Runge-Kutta 4, 14 stages version - Niegemann et al 2012
! Fastest RK4 scheme implemented, but less accurate then Carpenter RK4
! Very good performance for high N

CHARACTER(LEN=255),PARAMETER :: TimeDiscName = 'Niegemann RK4-14'
INTEGER,PARAMETER  :: nRKStages=14
REAL,PARAMETER  :: RK_a(2:14) = (/0.7188012108672410,&
                                   0.7785331173421570,&
                                   0.0053282796654044,&
                                   0.8552979934029281,&
                                   3.9564138245774565,&
                                   1.5780575380587385,&
                                   2.0837094552574054,&
                                   0.7483334182761610,&
                                   0.7032861106563359,&
                                  -0.0013917096117681,&
                                   0.0932075369637460,&
                                   0.9514200470875948,&
                                   7.1151571693922548 /)

REAL,PARAMETER  :: RK_b(1:14) = (/0.0367762454319673,&
                                   0.3136296607553959,&
                                   0.1531848691869027,&
                                   0.0030097086818182,&
                                   0.3326293790646110,&
                                   0.2440251405350864,&
                                   0.3718879239592277,&
                                   0.6204126221582444,&
                                   0.1524043173028741,&
                                   0.0760894927419266,&
                                   0.0077604214040978,&
                                   0.0024647284755382,&
                                   0.0780348340049386,&
                                   5.5059777270269628 /)

REAL,PARAMETER  :: RK_c(2:14) = (/0.0367762454319673,&
                                   0.1249685262725025,&
                                   0.2446177702277698,&
                                   0.2476149531070420,&
                                   0.2969311120382472,&
                                   0.3978149645802642,&
                                   0.5270854589440328,&
                                   0.6981269994175695,&
                                   0.8190890835352128,&
                                   0.8527059887098624,&
                                   0.8604711817462826,&
                                   0.8627060376969976,&
                                   0.8734213127600976 /)

! for degree N=1..15
REAL,PARAMETER  :: CFLScaleAlpha(1:15) = &
#if (PP_NodeType==1)
(/ 6.9716, 5.7724, 5.1863, 4.7880, 4.4741, 4.2120, 3.9836, 3.7811, 3.5617, 3.3705, 3.1995, 3.0488, 2.9137, 2.7900, 2.6786 /)
#elif (PP_NodeType==2)
(/ 14.7882, 9.5906, 7.9447, 7.0965, 6.5486, 6.1436, 5.8185, 5.5440, 5.3055, 5.0940, 4.9028, 4.7295, 4.5697, 4.4235, 4.2885 /)
#endif /*PP_NodeType*/
#endif
!===================================================================================================================================
END MODULE MOD_TimeDisc_Vars
