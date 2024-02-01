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
REAL             :: TEnd                              !> simulation end time
REAL             :: TAnalyze                          !> time for next analyze
REAL             :: Time                              !> Simulation Time
REAL             :: dt                                !> simulation time step
REAL             :: dtWeight                          !> part of original dt that is currently used as dt (output, rk, ...)
REAL             :: CFLScale                          !> cfl scale
REAL             :: CFLtoOne                          !> scaling factor to scale CFL to one
REAL             :: sdtCFLOne                         !> inverse of dt of CFLOne
!REAL             :: eps_LinearSolver,eps2_LinearSolver,epsTilde_LinearSolver
REAL             :: RKdtFrac,RKdtFracTotal
!INTEGER          :: maxIter_LinearSolver
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
#if (defined(IMPA) || defined(ROS) || (PP_TimeDiscMethod==509))
REAL             :: dt_old
#endif /*defined(IMPA) || defined(ROS) || (PP_TimeDiscMethod==509)*/
#if (PP_TimeDiscMethod==100)
INTEGER,PARAMETER:: nRKStages=1
#endif
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==2)|| (PP_TimeDiscMethod==6)
REAL,ALLOCATABLE :: Ut_temp(   :,:,:,:,:)             !> temporal variable for Ut
REAL,ALLOCATABLE :: U2t_temp(  :,:,:,:,:)             !> temporal variable for U2t
#ifdef PP_POIS
REAL,ALLOCATABLE :: Phit_temp( :,:,:,:,:)
#endif /*PP_POIS*/
#endif

!-----------------------------------------------------------------------------------------------------------------------------------
! TIME INTEGRATION: RUNGE_KUTTA COEFFICIENTS AND STABILITY NUMBERS
!-----------------------------------------------------------------------------------------------------------------------------------
#if (PP_TimeDiscMethod==1)||(PP_TimeDiscMethod==501)|| (PP_TimeDiscMethod==441)
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

#if (PP_TimeDiscMethod==6)||(PP_TimeDiscMethod==506)||(PP_TimeDiscMethod==443)
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

#if (PP_TimeDiscMethod==120)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! IMEX Heun/Crank-Nicolson1-2-2, implicit 1 stage, explicit 2 stages, common order 2  (A-stable, but not L-stable)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
INTEGER,PARAMETER :: nRKStages = 2
REAL,PARAMETER    :: ERK_a(2:nRKStages,1:nRKStages-1)= 1.0 ! RESHAPE( (/ 1.0, 0.0 /) , (/nRKStages-1,nRKStages/), ORDER=(/2,1/))
REAL,PARAMETER    :: RK_b(1:nRKStages) = (/ 0.5, 0.5 /)
REAL,PARAMETER    :: RK_c(2:nRKStages) = 1.0
REAL,PARAMETER    :: ESDIRK_a(2:nRKStages,1:nRKStages)= RESHAPE( (/ 0.5, 0.5 /) , (/nRKStages-1,nRKStages/), ORDER=(/2,1/))
! embedded scheme is euler
#endif
#if (PP_TimeDiscMethod==121)
! for stage counter
INTEGER,PARAMETER :: nRKStages = 4
! explicit scheme
! aij
REAL,PARAMETER  :: ERK3_a21=      1767732205903.0 /  2027836641118.0
REAL,PARAMETER  :: ERK3_a31=      5535828885825.0 / 10492691773637.0
REAL,PARAMETER  :: ERK3_a32=       788022342437.0 / 10882634858940.0
REAL,PARAMETER  :: ERK3_a41=      6485989280629.0 / 16251701735622.0
REAL,PARAMETER  :: ERK3_a42=     -4246266847089.0 /  9704473918619.0
REAL,PARAMETER  :: ERK3_a43=     10755448449292.0 / 10357097424841.0
REAL,PARAMETER  :: ERK3_a2(1:3) = (/ERK3_a21,      0.,      0./)
REAL,PARAMETER  :: ERK3_a3(1:3) = (/ERK3_a31,ERK3_a32,      0./)
REAL,PARAMETER  :: ERK3_a4(1:3) = (/ERK3_a41,ERK3_a42,ERK3_a43/)
REAL,PARAMETER  :: ERK_a(2:4,1:3) = RESHAPE( (/ERK3_a2(:),ERK3_a3(:),ERK3_a4(:)/),(/3,3/),ORDER =(/2,1/))
! both methods
! b i
REAL,PARAMETER  :: RK3_b1=  1471266399579.0 /   7840856788654.0
REAL,PARAMETER  :: RK3_b2= -4482444167858.0 /   7529755066697.0
REAL,PARAMETER  :: RK3_b3= 11266239266428.0 /  11593286722821.0
REAL,PARAMETER  :: RK3_b4=  1767732205903.0 /   4055673282236.0
REAL,PARAMETER  :: RK_b(1:4) = (/RK3_b1,RK3_b2,RK3_b3,RK3_b4/)
! bdach i maybe wrong
REAL,PARAMETER  :: RK3_bd1=   2756255671327.0 /  12835298489170.0
REAL,PARAMETER  :: RK3_bd2= -10771552573575.0 /   7529755066697.0
REAL,PARAMETER  :: RK3_bd3=   9247589265047.0 /  10645013368117.0
REAL,PARAMETER  :: RK3_bd4=   2193209047091.0 /   5459859503100.0
REAL,PARAMETER  :: RK_bd(1:4) = (/RK3_bd1,RK3_bd2,RK3_bd3,RK3_bd4/)
! c
REAL,PARAMETER  :: RK3_c2=   1767732205903.0  / 2027836641118.0
REAL,PARAMETER  :: RK3_c3=               3.0  / 5.0
REAL,PARAMETER  :: RK3_c4=               1.0
REAL,PARAMETER  :: RK_c(2:4) = (/RK3_c2,RK3_c3,RK3_c4/)
! implicit scheme
! aij
REAL,PARAMETER  :: ESDIRK3_a21=   1767732205903.0 /   4055673282236.0
REAL,PARAMETER  :: ESDIRK3_a22=   1767732205903.0 /   4055673282236.0
REAL,PARAMETER  :: ESDIRK3_a31=   2746238789719.0 /  10658868560708.0
REAL,PARAMETER  :: ESDIRK3_a32=   -640167445237.0 /   6845629431997.0
REAL,PARAMETER  :: ESDIRK3_a33=   1767732205903.0 /   4055673282236.0
REAL,PARAMETER  :: ESDIRK3_a41=   1471266399579.0 /   7840856788654.0
REAL,PARAMETER  :: ESDIRK3_a42=  -4482444167858.0 /   7529755066697.0
REAL,PARAMETER  :: ESDIRK3_a43=  11266239266428.0 /  11593286722821.0
REAL,PARAMETER  :: ESDIRK3_a44=   1767732205903.0 /   4055673282236.0
REAL,PARAMETER  :: ESDIRK3_a2(1:4) = (/ESDIRK3_a21,ESDIRK3_a22,         0.,         0./)
REAL,PARAMETER  :: ESDIRK3_a3(1:4) = (/ESDIRK3_a31,ESDIRK3_a32,ESDIRK3_a33,         0./)
REAL,PARAMETER  :: ESDIRK3_a4(1:4) = (/ESDIRK3_a41,ESDIRK3_a42,ESDIRK3_a43,ESDIRK3_a44/)
REAL,PARAMETER  :: ESDIRK_a(2:4,1:4)=RESHAPE((/ESDIRK3_a2,ESDIRK3_a3,ESDIRK3_a4/),(/3,4/),ORDER =(/2,1/))
! dense output
! for interpolation and extrapolation
! not implemented
REAL,PARAMETER :: RK3_bs11 =   4655552711362.0 /   22874653954995.0
REAL,PARAMETER :: RK3_bs12 =   -215264564351.0 /   13552729205753.0
REAL,PARAMETER :: RK3_bs1(1:2) = (/ RK3_bs11, RK3_bs12 /)
REAL,PARAMETER :: RK3_bs21 = -18682724506714.0 /    9892148508045.0
REAL,PARAMETER :: RK3_bs22 =  17890216137069.0 /   13817060693119.0
REAL,PARAMETER :: RK3_bs2(1:2) = (/ RK3_bs21, RK3_bs22 /)
REAL,PARAMETER :: RK3_bs31 =  34259539580243.0 /   13192909600954.0
REAL,PARAMETER :: RK3_bs32 = -28141676662227.0 /   17317692491321.0
REAL,PARAMETER :: RK3_bs3(1:2) = (/ RK3_bs31, RK3_bs32 /)
REAL,PARAMETER :: RK3_bs41 =    584795268549.0 /    6622622206610.0
REAL,PARAMETER :: RK3_bs42 =   2508943948391.0 /    7218656332882.0
REAL,PARAMETER :: RK3_bs4(1:2) = (/ RK3_bs41, RK3_bs42 /)
REAL,PARAMETER :: RK_bs(1:4,1:2)=RESHAPE((/RK3_bs1,RK3_bs2,RK3_bs3,RK3_bs4/),(/4,2/),ORDER =(/2,1/))
#endif
#if (PP_TimeDiscMethod==122)
! for stage counter
INTEGER,PARAMETER :: nRKStages = 6
! explicit scheme
! aij
REAL,PARAMETER  :: ERK4_a21= 1./2.
REAL,PARAMETER  :: ERK4_a31=           13861.0/62500.0
REAL,PARAMETER  :: ERK4_a32=            6889.0/62500.0
REAL,PARAMETER  :: ERK4_a41=  -116923316275.0 / 2393684061468.0
REAL,PARAMETER  :: ERK4_a42= -2731218467317.0 /15368042101831.0
REAL,PARAMETER  :: ERK4_a43=  9408046702089.0 /11113171139209.0
REAL,PARAMETER  :: ERK4_a51=  -451086348788.0 / 2902428689909.0
REAL,PARAMETER  :: ERK4_a52= -2682348792572.0 / 7519795681897.0
REAL,PARAMETER  :: ERK4_a53= 12662868775082.0 /11960479115383.0
REAL,PARAMETER  :: ERK4_a54=  3355817975965.0 /11060851509271.0
REAL,PARAMETER  :: ERK4_a61=   647845179188.0 / 3216320057751.0
REAL,PARAMETER  :: ERK4_a62=    73281519250.0 / 8382639484533.0
REAL,PARAMETER  :: ERK4_a63=   552539513391.0 / 3454668386233.0
REAL,PARAMETER  :: ERK4_a64=  3354512671639.0 / 8306763924573.0
REAL,PARAMETER  :: ERK4_a65=           4040.0 /         17871.0
REAL,PARAMETER  :: ERK4_a2(1:5) = (/ERK4_a21,      0.,      0.,      0.,      0./)
REAL,PARAMETER  :: ERK4_a3(1:5) = (/ERK4_a31,ERK4_a32,      0.,      0.,      0./)
REAL,PARAMETER  :: ERK4_a4(1:5) = (/ERK4_a41,ERK4_a42,ERK4_a43,      0.,      0./)
REAL,PARAMETER  :: ERK4_a5(1:5) = (/ERK4_a51,ERK4_a52,ERK4_a53,ERK4_a54,      0./)
REAL,PARAMETER  :: ERK4_a6(1:5) = (/ERK4_a61,ERK4_a62,ERK4_a63,ERK4_a64,ERK4_a65/)
REAL,PARAMETER  :: ERK_a(2:6,1:5) = RESHAPE( (/ERK4_a2(:),ERK4_a3(:),ERK4_a4(:),ERK4_a5(:),ERK4_a6(:)/),(/5,5/),ORDER =(/2,1/))
! both methods
! b i
REAL,PARAMETER  :: RK4_b1=  82889.0 / 524892.0
REAL,PARAMETER  :: RK4_b2=      0.0
REAL,PARAMETER  :: RK4_b3=  15625.0 /  83664.0
REAL,PARAMETER  :: RK4_b4=  69875.0 / 102672.0
REAL,PARAMETER  :: RK4_b5=  -2260.0 /   8211.0
REAL,PARAMETER  :: RK4_b6=      1.0 /      4.0
REAL,PARAMETER  :: RK_b(1:6) = (/RK4_b1,RK4_b2,RK4_b3,RK4_b4,RK4_b5,RK4_b6/)
! bdach i
REAL,PARAMETER  :: RK4_bd1= 4586570599.0  /  29645900160.0
REAL,PARAMETER  :: RK4_bd2=          0.0
REAL,PARAMETER  :: RK4_bd3=  178811875.0  /    945068544.0
REAL,PARAMETER  :: RK4_bd4=  814220225.0  /   1159782912.0
REAL,PARAMETER  :: RK4_bd5=   -3700637.0  /     11593932.0
REAL,PARAMETER  :: RK4_bd6=      61727.0  /       225920.0
REAL,PARAMETER  :: RK_bd(1:6) = (/RK4_bd1,RK4_bd2,RK4_bd3,RK4_bd4,RK4_bd5,RK4_bd6/)
! c
REAL,PARAMETER  :: RK4_c2=     1.0    /   2.0
REAL,PARAMETER  :: RK4_c3=    83.0    / 250.0
REAL,PARAMETER  :: RK4_c4=    31.0    /  50.0
REAL,PARAMETER  :: RK4_c5=    17.0    /  20.0
REAL,PARAMETER  :: RK4_c6=     1.0
REAL,PARAMETER  :: RK_c(2:6) = (/RK4_c2,RK4_c3,RK4_c4,RK4_c5,RK4_c6/)
! implicit scheme
! aij
REAL,PARAMETER  :: ESDIRK4_a21=              1.0 /            4.0
REAL,PARAMETER  :: ESDIRK4_a22=              1.0 /            4.0
REAL,PARAMETER  :: ESDIRK4_a31=           8611.0 /        62500.0
REAL,PARAMETER  :: ESDIRK4_a32=          -1743.0 /        31250.0
REAL,PARAMETER  :: ESDIRK4_a33=              1.0 /            4.0
REAL,PARAMETER  :: ESDIRK4_a41=        5012029.0 /     34652500.0
REAL,PARAMETER  :: ESDIRK4_a42=        -654441.0 /      2922500.0
REAL,PARAMETER  :: ESDIRK4_a43=         174375.0 /       388108.0
REAL,PARAMETER  :: ESDIRK4_a44=              1.0 /            4.0
REAL,PARAMETER  :: ESDIRK4_a51=    15267082809.0 / 155376265600.0
REAL,PARAMETER  :: ESDIRK4_a52=      -71443401.0 /    120774400.0
REAL,PARAMETER  :: ESDIRK4_a53=      730878875.0 /    902184768.0
REAL,PARAMETER  :: ESDIRK4_a54=        2285395.0 /      8070912.0
REAL,PARAMETER  :: ESDIRK4_a55=              1.0 /            4.0
REAL,PARAMETER  :: ESDIRK4_a61=          82889.0 /       524892.0
REAL,PARAMETER  :: ESDIRK4_a62=              0.0
REAL,PARAMETER  :: ESDIRK4_a63=          15625.0 /        83664.0
REAL,PARAMETER  :: ESDIRK4_a64=          69875.0 /       102672.0
REAL,PARAMETER  :: ESDIRK4_a65=          -2260.0 /         8211.0
REAL,PARAMETER  :: ESDIRK4_a66=              1.0 /            4.0
REAL,PARAMETER  :: ESDIRK4_a2(1:6) = (/ESDIRK4_a21,ESDIRK4_a22,         0.,         0.,         0.,         0./)
REAL,PARAMETER  :: ESDIRK4_a3(1:6) = (/ESDIRK4_a31,ESDIRK4_a32,ESDIRK4_a33,         0.,         0.,         0./)
REAL,PARAMETER  :: ESDIRK4_a4(1:6) = (/ESDIRK4_a41,ESDIRK4_a42,ESDIRK4_a43,ESDIRK4_a44,         0.,         0./)
REAL,PARAMETER  :: ESDIRK4_a5(1:6) = (/ESDIRK4_a51,ESDIRK4_a52,ESDIRK4_a53,ESDIRK4_a54,ESDIRK4_a55,         0./)
REAL,PARAMETER  :: ESDIRK4_a6(1:6) = (/ESDIRK4_a61,ESDIRK4_a62,ESDIRK4_a63,ESDIRK4_a64,ESDIRK4_a65,ESDIRK4_a66/)
REAL,PARAMETER  :: ESDIRK_a(2:6,1:6)=RESHAPE((/ESDIRK4_a2,ESDIRK4_a3,ESDIRK4_a4,ESDIRK4_a5,ESDIRK4_a6/),(/5,6/),ORDER =(/2,1/))
! dense output
! for interpolation and extrapolation
! third order implemented
REAL,PARAMETER :: RK4_bsO311 =   6943876665148.0 /  7220017795957.0
REAL,PARAMETER :: RK4_bsO312 =       -54480133.0 /       30881146.0
REAL,PARAMETER :: RK4_bsO313 =   6818779379841.0 /  7100303317025.0
REAL,PARAMETER :: RK4_bsO31(1:3) = (/ RK4_bsO311, RK4_bsO312, RK4_bsO313  /)
REAL,PARAMETER :: RK4_bsO321 =               0.0
REAL,PARAMETER :: RK4_bsO322 =               0.0
REAL,PARAMETER :: RK4_bsO323 =               0.0
REAL,PARAMETER :: RK4_bsO32(1:3) = (/ RK4_bsO321, RK4_bsO322, RK4_bsO323  /)
REAL,PARAMETER :: RK4_bsO331 =   7640104374378.0 /  9702883013639.0
REAL,PARAMETER :: RK4_bsO332 =       -11436875.0 /       14766696.0
REAL,PARAMETER :: RK4_bsO333 =   2173542590792.0 / 12501825683035.0
REAL,PARAMETER :: RK4_bsO33(1:3) = (/ RK4_bsO331, RK4_bsO332, RK4_bsO333  /)
REAL,PARAMETER :: RK4_bsO341 = -20649996744609.0 /  7521556579894.0
REAL,PARAMETER :: RK4_bsO342 =       174696575.0 /       18121608.0
REAL,PARAMETER :: RK4_bsO343 = -31592104683404.0 /  5083833661969.0
REAL,PARAMETER :: RK4_bsO34(1:3) = (/ RK4_bsO341, RK4_bsO342, RK4_bsO343  /)
REAL,PARAMETER :: RK4_bsO351 =   8854892464581.0 / 2390941311638.0
REAL,PARAMETER :: RK4_bsO352 =       -12120380.0 /        966161.0
REAL,PARAMETER :: RK4_bsO353 =  61146701046299.0 / 7138195549469.0
REAL,PARAMETER :: RK4_bsO35(1:3) = (/ RK4_bsO351, RK4_bsO352, RK4_bsO353  /)
REAL,PARAMETER :: RK4_bsO361 = -11397109935349.0 / 6675773540249.0
REAL,PARAMETER :: RK4_bsO362 =            3843.0 /           706.0
REAL,PARAMETER :: RK4_bsO363 = -17219254887155.0 / 4939391667607.0
REAL,PARAMETER :: RK4_bsO36(1:3) = (/ RK4_bsO361, RK4_bsO362, RK4_bsO363  /)
REAL,PARAMETER :: RK_bsO3(1:6,1:3)=RESHAPE((/RK4_bsO31,RK4_bsO32,RK4_bsO33,RK4_bsO34,RK4_bsO35,RK4_bsO36/),(/6,3/),ORDER =(/2,1/))
!! second order
REAL,PARAMETER :: RK4_bs11 =   5701579834848.0 /  6164663940925.0
REAL,PARAMETER :: RK4_bs12 =  -7364557999481.0 /   9602213853517.0
REAL,PARAMETER :: RK4_bs1(1:2) = (/ RK4_bs11, RK4_bs12 /)
REAL,PARAMETER :: RK4_bs21 =               0.0
REAL,PARAMETER :: RK4_bs22 =               0.0
REAL,PARAMETER :: RK4_bs2(1:2) = (/ RK4_bs21, RK4_bs22 /)
REAL,PARAMETER :: RK4_bs31 =  13131138058924.0 / 17779730471019.0
REAL,PARAMETER :: RK4_bs32 =  -6355522249597.0 / 11518083130066.0
REAL,PARAMETER :: RK4_bs3(1:2) = (/ RK4_bs31, RK4_bs32 /)
REAL,PARAMETER :: RK4_bs41 = -28096677048929.0 / 11161768239540.0
REAL,PARAMETER :: RK4_bs42 =  29755736407445.0 /  9305094404071.0
REAL,PARAMETER :: RK4_bs4(1:2) = (/ RK4_bs41, RK4_bs42 /)
REAL,PARAMETER :: RK4_bs51 =  42062433452849.0 / 11720557422164.0
REAL,PARAMETER :: RK4_bs52 = -38886896333129.0 / 10063858340160.0
REAL,PARAMETER :: RK4_bs5(1:2) = (/ RK4_bs51, RK4_bs52 /)
REAL,PARAMETER :: RK4_bs61 = -25841894007917.0 / 14894670528776.0
REAL,PARAMETER :: RK4_bs62 =  22142945955077.0 / 11155272088250.0
REAL,PARAMETER :: RK4_bs6(1:2) = (/ RK4_bs61, RK4_bs62 /)
REAL,PARAMETER :: RK_bs(1:6,1:2)=RESHAPE((/RK4_bs1,RK4_bs2,RK4_bs3,RK4_bs4,RK4_bs5,RK4_bs6/),(/6,2/),ORDER =(/2,1/))
#endif
#if (PP_TimeDiscMethod==123)
INTEGER,PARAMETER :: nRKStages = 2
! Ascher 1997 2 - 3 - 3
REAL,PARAMETER  :: RK3_gamma = (3.0 + SQRT(3.0))/6.0
! only implicit
! c
REAL,PARAMETER  :: RK3_c1=     RK3_gamma
REAL,PARAMETER  :: RK3_c2=     1.0 - RK3_gamma
REAL,PARAMETER  :: RK_c(1:2) = (/RK3_c1,RK3_c2/)
! implicit scheme - DIRK
! aij
REAL,PARAMETER  :: DIRK3_a11=   RK3_gamma
REAL,PARAMETER  :: DIRK3_a21=   1.0 - 2*RK3_gamma
REAL,PARAMETER  :: DIRK3_a22=   RK3_gamma
REAL,PARAMETER  :: DIRK3_a1(1:2) = (/DIRK3_a11,       0./)
REAL,PARAMETER  :: DIRK3_a2(1:2) = (/DIRK3_a21,DIRK3_a22/)
REAL,PARAMETER ::  ESDIRK(1:2,1:2)=RESHAPE((/DIRK3_a1,DIRK3_a2/),(/2,2/),ORDER =(/2,1/))
! bi
REAL,PARAMETER  :: RK3_b1=   0.5
REAL,PARAMETER  :: RK3_b2=   0.5
REAL,PARAMETER  :: RK3_b(1:2) = (/RK3_b1,RK3_b2/)
#endif
#ifdef IMPA
! || (PP_TimeDiscMethod==131)
REAL               :: RK_inc(2:nRKStages), RK_inflow(2:nRKStages),RK_fillSF
#endif /*IMPA*/
#if (PP_TimeDiscMethod==130)
! coefficients of Ianelli-Baker RO2-2
! Bassi-Paper
INTEGER,PARAMETER :: nRKStages = 2
REAL,PARAMETER  :: RK_gamma=0.5*(2.-sqrt(2.))
REAL,PARAMETER  :: RK_a21=  4.-8.*RK_gamma
REAL,PARAMETER  :: RK_a2(1:2) = (/RK_a21,      0./)
REAL,PARAMETER  :: RK_a(2:2,1:2) = RESHAPE( (/RK_a2(:)/),(/1,2/),ORDER =(/2,1/))
REAL,PARAMETER  :: RK_g21=  0.
REAL,PARAMETER  :: RK_g2(1:2) = (/RK_g21,      0./)
REAL,PARAMETER  :: RK_g(2:2,1:2) = RESHAPE( (/RK_g2(:)/),(/1,2/),ORDER =(/2,1/))
REAL,PARAMETER  :: RK_c2 = SUM(RK_a2)*RK_gamma
REAL,PARAMETER  :: RK_c(2:nRKStages) = (/RK_c2/)
REAL,PARAMETER  :: RK_b1= (1-1./(8.*RK_gamma))/RK_gamma
REAL,PARAMETER  :: RK_b2= 1./(8.*RK_gamma**2)
REAL,PARAMETER  :: RK_b(1:nRKStages) = (/RK_b1,RK_b2/)
#endif
#if (PP_TimeDiscMethod==131)
! coefficients of Lang-Verwer ROS3P (RO3-3)
! Bassi-Paper
INTEGER,PARAMETER :: nRKStages = 3
REAL,PARAMETER  :: RK_gamma=0.5  + sqrt(3.)/6.
REAL,PARAMETER  :: RK_a21=  1./RK_gamma ! gamma^-1
REAL,PARAMETER  :: RK_a2(1:nRKStages) = (/RK_a21,0., 0./)
REAL,PARAMETER  :: RK_a31=  RK_a21
REAL,PARAMETER  :: RK_a32=  0.
REAL,PARAMETER  :: RK_a3(1:nRKStages) = (/RK_a31,RK_a32, 0./)
REAL,PARAMETER  :: RK_a(2:nRKStages,1:nRKStages) = RESHAPE( (/RK_a2(:), RK_a3(:) /),(/nRKStages-1,nRKStages/),ORDER =(/2,1/))
REAL,PARAMETER  :: RK_g21=  -RK_a21*RK_a21 ! -gamma^-2
REAL,PARAMETER  :: RK_g2(1:nRKStages) = (/RK_g21,0.,0./)
REAL,PARAMETER  :: RK_g31= -RK_a21 *(1.+ RK_a21*(2.-1./(2.*RK_gamma))) ! -gamma^-1(1+gamma^-1*(2-1/(2gamma)))
REAL,PARAMETER  :: RK_g32=  -RK_a21*(2.-1./(2.*RK_gamma)) ! -gamma^-1 * 2-1/(2*gamma)
REAL,PARAMETER  :: RK_g3(1:nRKStages) = (/RK_g31,RK_g32,0./)
REAL,PARAMETER  :: RK_g(2:nRKStages,1:nRKStages) = RESHAPE( (/RK_g2(:), RK_g3(:) /),(/nRKStages-1,nRKStages/),ORDER =(/2,1/))
REAL,PARAMETER  :: RK_c2 = SUM(RK_a2)*RK_gamma
REAL,PARAMETER  :: RK_c3 = SUM(RK_a3)*RK_gamma
REAL,PARAMETER  :: RK_c(2:nRKStages) = (/RK_c2,RK_c3/)
REAL,PARAMETER  :: RK_b1= RK_a21 *(1. + RK_a21*(2./3. - 1./(6.*RK_gamma)))
REAL,PARAMETER  :: RK_b2= RK_a21 * (2./3. - 1./(6.*RK_gamma))
REAL,PARAMETER  :: RK_b3= 1./(3.*RK_gamma)
REAL,PARAMETER  :: RK_b(1:nRKStages) = (/RK_b1,RK_b2,RK_b3/)
#endif
#if (PP_TimeDiscMethod==132)
! Shampne ROS4 (RO4-4)
INTEGER,PARAMETER :: nRKStages = 4
! RK_a(i,j)
REAL,PARAMETER  :: RK_a21=    2.
REAL,PARAMETER  :: RK_a2(1:nRKStages) = (/RK_a21,0., 0., 0./)
REAL,PARAMETER  :: RK_a31=  48./25.
REAL,PARAMETER  :: RK_a32=   6./25.
REAL,PARAMETER  :: RK_a3(1:nRKStages) = (/RK_a31,RK_a32, 0., 0./)
REAL,PARAMETER  :: RK_a41=  48./25.
REAL,PARAMETER  :: RK_a42=   6./25.
REAL,PARAMETER  :: RK_a43=   0.
REAL,PARAMETER  :: RK_a4(1:nRKStages) = (/RK_a41,RK_a42, RK_a43, 0./)
REAL,PARAMETER  :: RK_a(2:nRKStages,1:nRKStages) = RESHAPE( (/RK_a2(:), RK_a3(:), RK_a4(:)/),(/nRKStages-1,nRKStages/),ORDER =(/2,1/))
! RK_c(i,j)
REAL,PARAMETER  :: RK_gamma=0.5
REAL,PARAMETER  :: RK_g21             =  -8.
REAL,PARAMETER  :: RK_g2(1:nRKStages) = (/RK_g21, 0.,0.,0./)
REAL,PARAMETER  :: RK_g31             =  372./25
REAL,PARAMETER  :: RK_g32             =  12./5
REAL,PARAMETER  :: RK_g3(1:nRKStages) = (/RK_g31, RK_g32,0.,0./)
REAL,PARAMETER  :: RK_g41             =  -112./125.
REAL,PARAMETER  :: RK_g42             =   -54./125.
REAL,PARAMETER  :: RK_g43             =    -2./5.
REAL,PARAMETER  :: RK_g4(1:nRKStages) = (/RK_g41, RK_g42,RK_g43,0./)
REAL,PARAMETER  :: RK_g(2:nRKStages,1:nRKStages) = RESHAPE( (/RK_g2(:), RK_g3(:), RK_g4(:)/),(/nRKStages-1,nRKStages/),ORDER =(/2,1/))
! RK_b(i)
REAL,PARAMETER  :: RK_c2 = SUM(RK_a2)*RK_gamma
REAL,PARAMETER  :: RK_c3 = SUM(RK_a3)*RK_gamma
REAL,PARAMETER  :: RK_c4 = SUM(RK_a4)*RK_gamma
REAL,PARAMETER  :: RK_c(2:nRKStages) = (/RK_c2,RK_c3,RK_c4/)
REAL,PARAMETER  :: RK_b1= 19./9.
REAL,PARAMETER  :: RK_b2= 0.5
REAL,PARAMETER  :: RK_b3= 25./108.
REAL,PARAMETER  :: RK_b4= 125./108.
REAL,PARAMETER  :: RK_b(1:nRKStages) = (/RK_b1,RK_b2,RK_b3,RK_b4/)
#endif
#if (PP_TimeDiscMethod==133)
! Steinebach RODASP  (RO4-6)
INTEGER,PARAMETER :: nRKStages = 6
! RK_a(i,j)
REAL,PARAMETER  :: RK_a21=  3.
REAL,PARAMETER  :: RK_a2(1:nRKStages) = (/RK_a21,0., 0., 0.,0.,0./)
REAL,PARAMETER  :: RK_a31=  1.831036793486759e0
REAL,PARAMETER  :: RK_a32=  4.955183967433795e-1
REAL,PARAMETER  :: RK_a3(1:nRKStages) = (/RK_a31,RK_a32, 0., 0.,0.,0./)
REAL,PARAMETER  :: RK_a41=  2.304376582692669e0
REAL,PARAMETER  :: RK_a42= -5.249275245743001e-2
REAL,PARAMETER  :: RK_a43= -1.176798761832782e0
REAL,PARAMETER  :: RK_a4(1:nRKStages) = (/RK_a41,RK_a42, RK_a43, 0.,0.,0./)
REAL,PARAMETER  :: RK_a51= -7.170454962423024e0
REAL,PARAMETER  :: RK_a52= -4.741636671481785e0
REAL,PARAMETER  :: RK_a53= -1.631002631330971e1
REAL,PARAMETER  :: RK_a54= -1.062004044111401e0
REAL,PARAMETER  :: RK_a5(1:nRKStages) = (/RK_a51,RK_a52, RK_a53, RK_a54,0.,0./)
REAL,PARAMETER  :: RK_a61= -7.170454962423024e0
REAL,PARAMETER  :: RK_a62= -4.741636671481785e0
REAL,PARAMETER  :: RK_a63= -1.631002631330971e1
REAL,PARAMETER  :: RK_a64= -1.062004044111401
REAL,PARAMETER  :: RK_a65=  1.0
REAL,PARAMETER  :: RK_a6(1:nRKStages) = (/RK_a61,RK_a62, RK_a63, RK_a64,RK_a65,0./)
REAL,PARAMETER::RK_a(2:nRKStages,1:nRKStages)=RESHAPE((/RK_a2,RK_a3,RK_a4,RK_a5,RK_a6/),(/nRKStages-1,nRKStages/),ORDER =(/2,1/))
! RK_c(i,j)
REAL,PARAMETER  :: RK_gamma=0.25
REAL,PARAMETER  :: RK_g21             = -1.2e1
REAL,PARAMETER  :: RK_g2(1:nRKStages) = (/RK_g21, 0.,0.,0.,0.,0./)
REAL,PARAMETER  :: RK_g31             = -8.791795173947035e0
REAL,PARAMETER  :: RK_g32             = -2.207865586973518e0
REAL,PARAMETER  :: RK_g3(1:nRKStages) = (/RK_g31, RK_g32,0.,0.,0.,0./)
REAL,PARAMETER  :: RK_g41             =  1.081793056857153e1
REAL,PARAMETER  :: RK_g42             =  6.780270611428266e0
REAL,PARAMETER  :: RK_g43             =  1.953485944642410e1
REAL,PARAMETER  :: RK_g4(1:nRKStages) = (/RK_g41, RK_g42,RK_g43,0.,0.,0./)
REAL,PARAMETER  :: RK_g51             =  3.419095006749676e1
REAL,PARAMETER  :: RK_g52             =  1.549671153725963e1
REAL,PARAMETER  :: RK_g53             =  5.474760875964130e1
REAL,PARAMETER  :: RK_g54             =  1.416005392148534e1
REAL,PARAMETER  :: RK_g5(1:nRKStages) = (/RK_g51, RK_g52,RK_g53,RK_g54,0.,0./)
REAL,PARAMETER  :: RK_g61             =  3.462605830930532e1
REAL,PARAMETER  :: RK_g62             =  1.530084976114473e1
REAL,PARAMETER  :: RK_g63             =  5.699955578662667e1
REAL,PARAMETER  :: RK_g64             =  1.840807009793095e1
REAL,PARAMETER  :: RK_g65             = -5.714285714285717e0
REAL,PARAMETER  :: RK_g6(1:nRKStages) = (/RK_g61, RK_g62,RK_g63,RK_g64,RK_g65,0./)
REAL,PARAMETER::RK_g(2:nRKStages,1:nRKStages)=RESHAPE((/RK_g2,RK_g3,RK_g4,RK_g5,RK_g6/),(/nRKStages-1,nRKStages/),ORDER =(/2,1/))
! RK_b(i)
REAL,PARAMETER  :: RK_c2 = SUM(RK_a2)*RK_gamma
REAL,PARAMETER  :: RK_c3 = SUM(RK_a3)*RK_gamma
REAL,PARAMETER  :: RK_c4 = SUM(RK_a4)*RK_gamma
REAL,PARAMETER  :: RK_c5 = SUM(RK_a5)*RK_gamma
REAL,PARAMETER  :: RK_c6 = SUM(RK_a6)*RK_gamma
REAL,PARAMETER  :: RK_c(2:nRKStages) = (/RK_c2,RK_c3,RK_c4,RK_c5,RK_c6/)
REAL,PARAMETER  :: RK_b1= -7.170454962423024e0
REAL,PARAMETER  :: RK_b2= -4.741636671481785e0
REAL,PARAMETER  :: RK_b3= -1.631002631330971e1
REAL,PARAMETER  :: RK_b4= -1.062004044111401e0
REAL,PARAMETER  :: RK_b5=  1.0
REAL,PARAMETER  :: RK_b6=  1.0
REAL,PARAMETER  :: RK_b(1:nRKStages) = (/RK_b1,RK_b2,RK_b3,RK_b4,RK_b5,RK_b6/)
#endif
#if (PP_TimeDiscMethod==134)
! Kaps and Wanner ROW64 ROS6 (RO6-6)
INTEGER,PARAMETER :: nRKStages = 6
! RK_a(i,j)
REAL,PARAMETER  :: RK_a21= 2.
REAL,PARAMETER  :: RK_a2(1:nRKStages) = (/RK_a21,0., 0., 0.,0.,0./)
REAL,PARAMETER  :: RK_a31=  1.751493065942685e0
REAL,PARAMETER  :: RK_a32= -1.454290536332865e-1
REAL,PARAMETER  :: RK_a3(1:nRKStages) = (/RK_a31,RK_a32, 0., 0.,0.,0./)
REAL,PARAMETER  :: RK_a41= -1.847093912231436e0
REAL,PARAMETER  :: RK_a42= -2.513756792158473e0
REAL,PARAMETER  :: RK_a43=  1.874707432337999e0
REAL,PARAMETER  :: RK_a4(1:nRKStages) = (/RK_a41,RK_a42, RK_a43, 0.,0.,0./)
REAL,PARAMETER  :: RK_a51=  1.059634783677141e1
REAL,PARAMETER  :: RK_a52=  1.974951525952609e0
REAL,PARAMETER  :: RK_a53= -1.905211286263863e0
REAL,PARAMETER  :: RK_a54= -3.575118228830491e0
REAL,PARAMETER  :: RK_a5(1:nRKStages) = (/RK_a51,RK_a52, RK_a53, RK_a54,0.,0./)
REAL,PARAMETER  :: RK_a61=  2.417642067883312e0
REAL,PARAMETER  :: RK_a62=  3.050984437044573e-1
REAL,PARAMETER  :: RK_a63= -2.346208879122501e-1
REAL,PARAMETER  :: RK_a64= -1.327038464607418e-1
REAL,PARAMETER  :: RK_a65=  3.912922779645768e-2
REAL,PARAMETER  :: RK_a6(1:nRKStages) = (/RK_a61,RK_a62, RK_a63, RK_a64,RK_a65,0./)
REAL,PARAMETER::RK_a(2:nRKStages,1:nRKStages)=RESHAPE((/RK_a2,RK_a3,RK_a4,RK_a5,RK_a6/),(/nRKStages-1,nRKStages/),ORDER =(/2,1/))
! RK_c(i,j)
REAL,PARAMETER  :: RK_gamma=3.3414236706805040e-1
REAL,PARAMETER  :: RK_g21             = -1.7450294925129950e1
REAL,PARAMETER  :: RK_g2(1:nRKStages) = (/RK_g21, 0.,0.,0.,0.,0./)
REAL,PARAMETER  :: RK_g31             = -1.2023599362278440e1
REAL,PARAMETER  :: RK_g32             =  1.3159101107427450e0
REAL,PARAMETER  :: RK_g3(1:nRKStages) = (/RK_g31, RK_g32,0.,0.,0.,0./)
REAL,PARAMETER  :: RK_g41             =  2.3112305971592720e1
REAL,PARAMETER  :: RK_g42             =  1.2978931295654450e1
REAL,PARAMETER  :: RK_g43             = -8.4453745945620380e0
REAL,PARAMETER  :: RK_g4(1:nRKStages) = (/RK_g41, RK_g42,RK_g43,0.,0.,0./)
REAL,PARAMETER  :: RK_g51             = -3.1472288913307130e0
REAL,PARAMETER  :: RK_g52             = -1.7613326229099650e0
REAL,PARAMETER  :: RK_g53             =  6.1152959340385850e0
REAL,PARAMETER  :: RK_g54             =  1.4993199504571120e1
REAL,PARAMETER  :: RK_g5(1:nRKStages) = (/RK_g51, RK_g52,RK_g53,RK_g54,0.,0./)
REAL,PARAMETER  :: RK_g61             = -2.0158409112628800e1
REAL,PARAMETER  :: RK_g62             = -1.6039237998001330e0
REAL,PARAMETER  :: RK_g63             =  1.1558700969202520e0
REAL,PARAMETER  :: RK_g64             =  6.3046398152920440e-1
REAL,PARAMETER  :: RK_g65             = -1.6025102156371740e-1
REAL,PARAMETER  :: RK_g6(1:nRKStages) = (/RK_g61, RK_g62,RK_g63,RK_g64,RK_g65,0./)
REAL,PARAMETER::RK_g(2:nRKStages,1:nRKStages)=RESHAPE((/RK_g2,RK_g3,RK_g4,RK_g5,RK_g6/),(/nRKStages-1,nRKStages/),ORDER =(/2,1/))
! RK_b(i)
REAL,PARAMETER  :: RK_c1             = 0.
REAL,PARAMETER  :: RK_c2             = 0.66828473413610087e0
REAL,PARAMETER  :: RK_c3             = 0.82000000000000000e0
REAL,PARAMETER  :: RK_c4             = 0.21963625075792513e0
REAL,PARAMETER  :: RK_c5             = 0.89999999999999999e0
REAL,PARAMETER  :: RK_c6             = 0.66585763293194957e0
REAL,PARAMETER  :: RK_c(2:nRKStages) = (/RK_c2,RK_c3,RK_c4,RK_c5,RK_c6/)
REAL,PARAMETER  :: RK_b1= 3.3993474526741650e1
REAL,PARAMETER  :: RK_b2=-2.0918298828473330e1
REAL,PARAMETER  :: RK_b3=-1.3756884774710810e1
REAL,PARAMETER  :: RK_b4=-1.1139259299300770e1
REAL,PARAMETER  :: RK_b5= 2.8734065276094680e0
REAL,PARAMETER  :: RK_b6= 3.8766099456208400e1
REAL,PARAMETER  :: RK_b(1:nRKStages) = (/RK_b1,RK_b2,RK_b3,RK_b4,RK_b5,RK_b6/)
#endif
#ifdef ROS
REAL            :: dt_inv
REAL            :: RK_inflow(2:nRKStages) ! required for boundary conditions
#endif /*ROSENBROCK RK*/
!===================================================================================================================================
END MODULE MOD_TimeDisc_Vars
