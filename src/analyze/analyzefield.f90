#include "boltzplatz.h"

MODULE MOD_AnalyzeField
!===================================================================================================================================
! Contains the Poynting Vector Integral part for the power analysis of the field vector
!===================================================================================================================================
USE MOD_Globals, ONLY:UNIT_stdout
USE MOD_PreProc
!===================================================================================================================================
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE GetPoyntingIntPlane
  MODULE PROCEDURE GetPoyntingIntPlane
END INTERFACE

INTERFACE FinalizePoyntingInt
  MODULE PROCEDURE FinalizePoyntingInt
END INTERFACE

#if (PP_nVar>=6)
INTERFACE CalcPoyntingIntegral
  MODULE PROCEDURE CalcPoyntingIntegral
END INTERFACE
#endif

INTERFACE CalcPotentialEnergy
  MODULE PROCEDURE CalcPotentialEnergy
END INTERFACE

PUBLIC:: GetPoyntingIntPlane,FinalizePoyntingInt,CalcPotentialEnergy
#if (PP_nVar>=6)
PUBLIC:: CalcPoyntingIntegral
#endif
#ifndef PARTICLES
PUBLIC:: AnalyzeField
#endif /*NOT PARTICLES*/
!===================================================================================================================================

CONTAINS

#ifndef PARTICLES
SUBROUTINE AnalyzeField(Time)
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Analyze_Vars           ,ONLY: DoAnalyze
USE MOD_Particle_Analyze_Vars  ,ONLY: CalcEpot, IsRestart
USE MOD_Restart_Vars           ,ONLY: DoRestart
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: Time
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL             :: isOpen, FileExists                                          !
CHARACTER(LEN=350)  :: outfile                                                      !
INTEGER             :: unit_index, OutputCounter
REAL                :: WEl, WMag
!===================================================================================================================================
IF ( DoRestart ) THEN
  isRestart = .true.
END IF
IF (DoAnalyze) THEN
!SWRITE(UNIT_StdOut,'(132("-"))')
!SWRITE(UNIT_stdOut,'(A)') ' PERFORMING PARTICLE ANALYZE...'
OutputCounter = 2
unit_index = 535
#ifdef MPI
 IF(MPIROOT)THEN
#endif    /* MPI */
    INQUIRE(UNIT   = unit_index , OPENED = isOpen)
    IF (.NOT.isOpen) THEN
      outfile = 'Database.csv'
      INQUIRE(file=TRIM(outfile),EXIST=FileExists)
      IF (isRestart .and. FileExists) THEN
         OPEN(unit_index,file=TRIM(outfile),position="APPEND",status="OLD")
         !CALL FLUSH (unit_index)
      ELSE
         OPEN(unit_index,file=TRIM(outfile))
         !CALL FLUSH (unit_index)
         !--- insert header
       
         WRITE(unit_index,'(A6,A5)',ADVANCE='NO') 'TIME', ' '
         IF (CalcEpot) THEN 
           WRITE(unit_index,'(A1)',ADVANCE='NO') ','
           WRITE(unit_index,'(I3.3,A11)',ADVANCE='NO') OutputCounter,'-W-El      '
             OutputCounter = OutputCounter + 1
           WRITE(unit_index,'(A1)',ADVANCE='NO') ','
           WRITE(unit_index,'(I3.3,A11)',ADVANCE='NO') OutputCounter,'-W-Mag    '
             OutputCounter = OutputCounter + 1
         END IF
         WRITE(unit_index,'(A14)') ' ' 
      END IF
    END IF
#ifdef MPI
 END IF
#endif    /* MPI */


!IF (CalcCharge.AND.(.NOT.ChargeCalcDone)) CALL CalcDepositedCharge()
IF(CalcEpot) CALL CalcPotentialEnergy(WEl,WMag)

#ifdef MPI
 IF(MPIROOT)THEN
#endif    /* MPI */
   WRITE(unit_index,104,ADVANCE='NO') Time
   IF (CalcEpot) THEN 
     WRITE(unit_index,'(A1)',ADVANCE='NO') ','
     WRITE(unit_index,104,ADVANCE='NO') WEl
     WRITE(unit_index,'(A1)',ADVANCE='NO') ','
     WRITE(unit_index,104,ADVANCE='NO') WMag
   END IF
   WRITE(unit_index,'(A1)') ' ' 
#ifdef MPI
 END IF
#endif    /* MPI */

104    FORMAT (e25.14)

!SWRITE(UNIT_stdOut,'(A)')' PARTCILE ANALYZE DONE!'
!SWRITE(UNIT_StdOut,'(132("-"))')
ELSE
!SWRITE(UNIT_stdOut,'(A)')' NO PARTCILE ANALYZE TO DO!'
!SWRITE(UNIT_StdOut,'(132("-"))')
END IF ! DoAnalyze

END SUBROUTINE AnalyzeField
#endif /*NOT PARTICLES*/

#if (PP_nVar>=6)
SUBROUTINE CalcPoyntingIntegral(t,doProlong)
!===================================================================================================================================
! Calculation of Poynting Integral with its own Prolong to face // check if Gauss-Labatto or Gaus Points is used is missing ... ups
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars             ,ONLY:isPoyntingIntSide,nElems, SurfElem, NormVec,whichPoyntingPlane
USE MOD_Mesh_Vars             ,ONLY:ElemToSide,PoyntingMainDir
USE MOD_Analyze_Vars          ,ONLY:nPoyntingIntPlanes, S!, STEM
USE MOD_Interpolation_Vars    ,ONLY:L_Minus,L_Plus,wGPSurf
USE MOD_DG_Vars               ,ONLY:U,U_master
USE MOD_Equation_Vars         ,ONLY:smu0
#ifdef MPI
  USE MOD_Globals
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)          :: t
LOGICAL,INTENT(IN),OPTIONAL :: doProlong
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: iElem, SideID,ilocSide,iPoyntingSide
INTEGER          :: p,q,l
REAL             :: Uface(PP_nVar,0:PP_N,0:PP_N)
REAL             :: SIP(0:PP_N,0:PP_N)
REAL             :: Sabs(nPoyntingIntPlanes), STEMabs(nPoyntingIntPlanes)
#ifdef MPI
REAL             :: SumSabs(nPoyntingIntPlanes)
#endif
LOGICAL          :: Prolong=.TRUE.
!REAL             :: sresvac
!===================================================================================================================================

IF(PRESENT(doProlong))THEN
  Prolong=doProlong
ELSE
  Prolong=.TRUE.
ENDIF  
! TEM coefficient
!sresvac = 1./sqrt(mu0/eps0)

S    = 0.
!STEM = 0.
Sabs = 0.
STEMabs = 0.

iPoyntingSide = 0 ! only if all poynting vectors are desired
DO iELEM = 1, nElems
  Do ilocSide = 1, 6
    IF(ElemToSide(E2S_FLIP,ilocSide,iElem)==0)THEN ! only master sides
      SideID=ElemToSide(E2S_SIDE_ID,ilocSide,iElem)
      IF(.NOT.isPoyntingIntSide(SideID)) CYCLE
      IF(Prolong)THEN
#if (PP_NodeType==1) /* for Gauss-points*/
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,0,p,q,iElem)*L_Minus(0)
              DO l=1,PP_N
                ! switch to right hand system
                Uface(:,q,p)=Uface(:,q,p)+U(:,l,p,q,iElem)*L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ETA_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,p,0,q,iElem)*L_Minus(0)
              DO l=1,PP_N
                Uface(:,p,q)=Uface(:,p,q)+U(:,p,l,q,iElem)*L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ZETA_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,p,q,0,iElem)*L_Minus(0)
              DO l=1,PP_N
                ! switch to right hand system
                Uface(:,q,p)=Uface(:,q,p)+U(:,p,q,l,iElem)*L_Minus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(XI_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,0,p,q,iElem)*L_Plus(0)
              DO l=1,PP_N
                Uface(:,p,q)=Uface(:,p,q)+U(:,l,p,q,iElem)*L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,PP_N-p,q)=U(:,p,0,q,iElem)*L_Plus(0)
              DO l=1,PP_N
                ! switch to right hand system
                Uface(:,PP_N-p,q)=Uface(:,PP_N-p,q)+U(:,p,l,q,iElem)*L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        CASE(ZETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,p,q,0,iElem)*L_Plus(0)
              DO l=1,PP_N
                Uface(:,p,q)=Uface(:,p,q)+U(:,p,q,l,iElem)*L_Plus(l)
              END DO ! l
            END DO ! p
          END DO ! q
        END SELECT
#else /* for Gauss-Lobatto-points*/
        SELECT CASE(ilocSide)
        CASE(XI_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,0,p,q,iElem)
            END DO ! p
          END DO ! q
        CASE(ETA_MINUS)
          Uface(:,:,:)=U(:,:,0,:,iElem)
        CASE(ZETA_MINUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,q,p)=U(:,p,q,0,iElem)
            END DO ! p
          END DO ! q
        CASE(XI_PLUS)
          Uface(:,:,:)=U(:,PP_N,:,:,iElem)
        CASE(ETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,PP_N-p,q)=U(:,p,PP_N,q,iElem)
            END DO ! p
          END DO ! q
        CASE(ZETA_PLUS)
          DO q=0,PP_N
            DO p=0,PP_N
              Uface(:,p,q)=U(:,p,q,PP_N,iElem)
            END DO ! p
          END DO ! q
        END SELECT
#endif
        ELSE ! no prolonge to face
          Uface=U_master(:,:,:,SideID)
        END IF ! Prolong
        ! calculate poynting vector
        iPoyntingSide = iPoyntingSide + 1
        CALL PoyntingVector(Uface(:,:,:),S(:,:,:,iPoyntingSide))
        IF ( NormVec(PoyntingMainDir,0,0,SideID) .GT. 0 ) THEN
          SIP(:,:) = S(1,:,:,iPoyntingSide) * NormVec(1,:,:,SideID) &
                   + S(2,:,:,iPoyntingSide) * NormVec(2,:,:,SideID) &
                   + S(3,:,:,iPoyntingSide) * NormVec(3,:,:,SideID)
        ELSE ! NormVec(PoyntingMainDir,:,:,iPoyningSide) < 0
          SIP(:,:) =-S(1,:,:,iPoyntingSide) * NormVec(1,:,:,SideID) &
                   - S(2,:,:,iPoyntingSide) * NormVec(2,:,:,SideID) &
                   - S(3,:,:,iPoyntingSide) * NormVec(3,:,:,SideID)
        END IF ! NormVec(PoyntingMainDir,:,:,iPoyntingSide)
        ! multiplied by surface element and  Gaus Points
        SIP(:,:) = SIP(:,:) * SurfElem(:,:,SideID) * wGPSurf(:,:)
        ! total flux through each plane
        Sabs(whichPoyntingPlane(SideID)) = Sabs(whichPoyntingPlane(SideID)) + smu0* SUM(SIP(:,:))
    END IF ! flip =0
  END DO ! iSides
END DO ! iElems

#ifdef MPI
  CALL MPI_REDUCE   (Sabs(:) , sumSabs(:) , nPoyntingIntPlanes , MPI_DOUBLE_PRECISION ,MPI_SUM, 0, MPI_COMM_WORLD,IERROR)
  Sabs(:) = sumSabs(:)
#endif /* MPI */

! output callling
CALL OutputPoyntingInt(t,Sabs(:)) 

END SUBROUTINE CalcPoyntingIntegral
#endif


#if (PP_nVar>=6)
SUBROUTINE PoyntingVector(Uface_in,Sloc)
!===================================================================================================================================
! Calculate the Poynting Vector on a certain face
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)       :: Uface_in(PP_nVar,0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)      :: Sloc(1:3,0:PP_N,0:PP_N)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: p,q
!===================================================================================================================================

! calculate the poynting vector at each node, additionally the abs of the poynting vector only based on E
DO p = 0,PP_N
  DO q = 0,PP_N
    Sloc(1,p,q)  =  Uface_in(2,p,q)*Uface_in(6,p,q) - Uface_in(3,p,q)*Uface_in(5,p,q) 
    Sloc(2,p,q)  = -Uface_in(1,p,q)*Uface_in(6,p,q) + Uface_in(3,p,q)*Uface_in(4,p,q) 
    Sloc(3,p,q)  =  Uface_in(1,p,q)*Uface_in(5,p,q) - Uface_in(2,p,q)*Uface_in(4,p,q) 
  END DO ! q - PP_N
END DO  ! p - PP_N

END SUBROUTINE PoyntingVector
#endif


#if (PP_nVar>=6)
SUBROUTINE OutputPoyntingInt(t,Sabs)
!===================================================================================================================================
! Output of PoyntingVector Integral to *csv vile
!===================================================================================================================================
! MODULES
USE MOD_Analyze_Vars          ,ONLY:nPoyntingIntPlanes,PosPoyntingInt
USE MOD_Restart_Vars          ,ONLY:DoRestart
#ifdef MPI
  USE MOD_Globals
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: t, Sabs(nPoyntingIntPlanes)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: unit_index_PI, iPlane
LOGICAL             :: isRestart, isOpen,FileExists
CHARACTER(LEN=64)   :: filename_PI
!===================================================================================================================================
isRestart=.FALSE.
IF (DoRestart) THEN
  isRestart=.TRUE.
END IF

filename_PI  = 'Power.csv'
unit_index_PI=273

#ifdef MPI
IF(MPIRoot)THEN
#endif    /* MPI */

INQUIRE(UNIT   = unit_index_PI , OPENED = isOpen)
IF (.NOT.isOpen) THEN
  INQUIRE(file=TRIM(filename_PI),EXIST=FileExists)
  IF (isRestart .and. FileExists) THEN
    OPEN(unit_index_PI,file=TRIM(filename_PI),position="APPEND",status="OLD")
  ELSE
    OPEN(unit_index_PI,file=TRIM(filename_PI))
    ! --- insert header
    WRITE(unit_index_PI,'(A6,A5)',ADVANCE='NO') 'TIME', ' '
    DO iPlane = 1, nPoyntingIntPlanes
      WRITE(unit_index_PI,'(A1)',ADVANCE='NO') ','
      WRITE(unit_index_PI,'(A14,I0.3,A1,E13.7,A1)',ADVANCE='NO') 'Plane-Pos-',iPlane,'(', PosPoyntingInt(iPlane),')'
    END DO              
    WRITE(unit_index_PI,'(A1)') ''
  END IF
END IF
! write data to file
WRITE(unit_index_PI,'(e25.14)',ADVANCE='NO') t
DO iPlane = 1, nPoyntingIntPlanes
  WRITE(unit_index_PI,'(A1)',ADVANCE='NO') ','
  WRITE(unit_index_PI,'(e25.14)',ADVANCE='NO') Sabs(iPlane)
END DO
WRITE(unit_index_PI,'(A1)') ''

#ifdef MPI
 END IF
#endif    /* MPI */

END SUBROUTINE OutputPoyntingInt
#endif

SUBROUTINE GetPoyntingIntPlane()
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars         ,ONLY:nPoyntingIntSides, isPoyntingIntSide,nSides,nElems,Face_xGP,whichPoyntingPlane
USE MOD_Mesh_Vars         ,ONLY:ElemToSide,normvec,PoyntingMainDir
USE MOD_Analyze_Vars      ,ONLY:PoyntingIntCoordErr,nPoyntingIntPlanes,PosPoyntingInt,PoyntingIntPlaneFactor , S, STEM
USE MOD_ReadInTools       ,ONLY:GETINT,GETREAL
#ifdef MPI
  USE MOD_Globals
  !USE MOD_part_MPI_Vars   ,ONLY:PMPIVAR
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem, iSide, iPlane, SideID
INTEGER,ALLOCATABLE :: nFaces(:)
REAL                :: diff
INTEGER             :: p,q
CHARACTER(LEN=32)   :: index_plane
INTEGER,ALLOCATABLE :: sumFaces(:)
INTEGER             :: sumAllfaces
INTEGER             :: PoyntingNormalDir1,PoyntingNormalDir2
!===================================================================================================================================

SWRITE(UNIT_stdOut,'(A)') ' GET PLANES TO CALCULATE POYNTING VECTOR INTEGRAL ...'

! first stuff
nPoyntingIntSides=0 
ALLOCATE(isPoyntingIntSide(1:nSides))
isPoyntingIntSide = .FALSE.

! know get number of planes and coordinates
nPoyntingIntPlanes   = GETINT('PoyntingVecInt-Planes','0')
PoyntingMainDir = GETINT('PoyntingMainDir','3') ! default "3" is z-direction 
SELECT CASE (PoyntingMainDir)
  CASE (1) ! poynting vector integral in x-direction
    PoyntingNormalDir1=2
    PoyntingNormalDir2=3
  CASE (2) ! poynting vector integral in y-direction
    PoyntingNormalDir1=1
    PoyntingNormalDir2=3
  CASE (3) ! poynting vector integral in z-direction
    PoyntingNormalDir1=1
    PoyntingNormalDir2=2
  CASE DEFAULT
    CALL abort(&
    __STAMP__&
    ,'Poynting vector itnegral currently only in x,y,z!')
END SELECT
ALLOCATE(PosPoyntingInt(nPoyntingIntPlanes))
ALLOCATE(PoyntingIntPlaneFactor(nPoyntingIntPlanes))
ALLOCATE(whichPoyntingPlane(nSides))
ALLOCATE(nFaces(nPoyntingIntPlanes))
whichPoyntingPlane = -1
nFaces(:) = 0

DO iPlane=1,nPoyntingIntPlanes
 WRITE(UNIT=index_plane,FMT='(I2.2)') iPlane 
  SELECT CASE (PoyntingMainDir)
    CASE (1)
      PosPoyntingInt(iPlane)= GETREAL('Plane-'//TRIM(index_plane)//'-x-coord','0.')
    CASE (2)
      PosPoyntingInt(iPlane)= GETREAL('Plane-'//TRIM(index_plane)//'-y-coord','0.')
    CASE (3)
      PosPoyntingInt(iPlane)= GETREAL('Plane-'//TRIM(index_plane)//'-z-coord','0.')
END SELECT
 PoyntingIntPlaneFactor= GETREAL('Plane-'//TRIM(index_plane)//'-factor','1.')
END DO
PoyntingIntCoordErr=GETREAL('Plane-Tolerance','1E-5')

! loop over all planes
DO iPlane = 1, nPoyntingIntPlanes
  ! loop over all elements
  DO iElem=1,nElems
    ! loop over all local sides
    DO iSide=1,6
      IF(ElemToSide(E2S_FLIP,iSide,iElem)==0)THEN ! only master sides
        SideID=ElemToSide(E2S_SIDE_ID,iSide,iElem)
        ! first search only planes with normal vector parallel to direction of "MainDir"
        IF((     NormVec(PoyntingNormalDir1,0,0,SideID)  < PoyntingIntCoordErr) .AND. &
           (     NormVec(PoyntingNormalDir2,0,0,SideID)  < PoyntingIntCoordErr) .AND. &
           ( ABS(NormVec(PoyntingMainDir   ,0,0,SideID)) > PoyntingIntCoordErr))THEN
        ! loop over all Points on Face
          DO q=0,PP_N
            DO p=0,PP_N
              diff = ABS(Face_xGP(PoyntingMainDir,p,q,SideID) - PosPoyntingInt(iPlane))
              IF (diff < PoyntingIntCoordErr) THEN
                IF (.NOT.isPoyntingIntSide(SideID)) THEN
                  nPoyntingIntSides = nPoyntingIntSides +1
                  whichPoyntingPlane(SideID) = iPlane
                  isPoyntingIntSide(SideID) = .TRUE.
                  nFaces(iPlane) = nFaces(iPlane) + 1
                END IF
              END IF ! diff < eps
            END DO !p
          END DO !q
        END IF ! n parallel gyrotron axis
      END IF ! flip = 0 master side
    END DO ! iSides
  END DO !iElem=1,nElems
END DO ! iPlanes

ALLOCATE(sumFaces(nPoyntingIntPlanes))
#ifdef MPI
sumFaces=0
sumAllFaces=0
  CALL MPI_REDUCE(nFaces , sumFaces , nPoyntingIntPlanes , MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, IERROR)
  !nFaces(:) = sumFaces(:)
  CALL MPI_REDUCE(nPoyntingIntSides , sumAllFaces , 1 , MPI_INTEGER, MPI_SUM,0, MPI_COMM_WORLD, IERROR)
  !nPoyntingIntSides = sumAllFaces
#else
sumFaces=nFaces
sumAllFaces=nPoyntingIntSides
#endif /* MPI */

DO iPlane= 1, nPoyntingIntPlanes
  SWRITE(UNIT_stdOut,'(A,I2,A,I10,A)') 'Processed plane no.: ',iPlane,'. Found ',sumFaces(iPlane),' surfaces.'
END DO
SWRITE(UNIT_stdOut,'(A,I10,A)') 'A total of',sumAllFaces, &
                        ' surfaces for the poynting vector integral calculation are found.'

ALLOCATE(S    (1:3,0:PP_N,0:PP_N,1:nPoyntingIntSides) , &
         STEM     (0:PP_N,0:PP_N,1:nPoyntingIntSides)  )

SWRITE(UNIT_stdOut,'(A)') ' ... POYNTING VECTOR INTEGRAL INITIALIZATION DONE.'  

END SUBROUTINE GetPoyntingIntPlane

SUBROUTINE FinalizePoyntingInt()
!===================================================================================================================================
! Finalize Poynting Integral
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars         ,ONLY:isPoyntingIntSide,whichPoyntingPlane
USE MOD_Analyze_Vars      ,ONLY:PosPoyntingInt,PoyntingIntPlaneFactor, S, STEM
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! DEALLOCATE ALL
SDEALLOCATE(isPoyntingIntSide)
SDEALLOCATE(PosPoyntingInt)
SDEALLOCATE(PoyntingIntPlaneFactor)
SDEALLOCATE(whichPoyntingPlane)
SDEALLOCATE(S)
SDEALLOCATE(STEM)

END SUBROUTINE FinalizePoyntingInt

SUBROUTINE CalcPotentialEnergy(WEl, WMag) 
!===================================================================================================================================
! Initializes variables necessary for analyse subroutines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,          ONLY : nElems, sJ
USE MOD_Interpolation_Vars, ONLY : wGP
USE MOD_Equation_Vars,      ONLY : smu0, eps0 
#ifndef PP_HDG
USE MOD_DG_Vars,            ONLY : U
USE MOD_Mesh_Vars,          ONLY : Elem_xGP
#endif /*PP_nVar=8*/        
#ifdef PP_HDG
#if PP_nVar==1
USE MOD_Equation_Vars,        ONLY:E
#elif PP_nVar==3
USE MOD_Equation_Vars,        ONLY:B
#else
USE MOD_Equation_Vars,        ONLY:B,E
#endif /*PP_nVar==1*/
#else
USE MOD_PML_Vars,           ONLY : xyzPhysicalMinMax,DoPML
#endif /*PP_HDG*/
#ifdef MPI
#endif /*MPI*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: WEl, WMag 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iElem
INTEGER           :: i,j,k
REAL              :: J_N(1,0:PP_N,0:PP_N,0:PP_N)
REAL              :: WEl_tmp, WMag_tmp, E_abs
#ifndef PP_HDG
REAL              :: B_abs 
#endif
#ifdef MPI
REAL              :: RD
#endif
!===================================================================================================================================

Wel=0.
WMag=0.

#ifndef PP_HDG
IF(DoPML)THEN
  DO iElem=1,nElems
    !--- Calculate and save volume of element iElem
    WEl_tmp=0. 
    WMag_tmp=0. 
    J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  ! in electromagnetische felder by henke 2011 - springer
  ! WMag = 1/(2mu) * int_V B^2 dV 
      E_abs = U(1,i,j,k,iElem)*U(1,i,j,k,iElem) &
            + U(2,i,j,k,iElem)*U(2,i,j,k,iElem) &
            + U(3,i,j,k,iElem)*U(3,i,j,k,iElem)
#if (PP_nVar==8)
      B_abs = U(4,i,j,k,iElem)*U(4,i,j,k,iElem) &
            + U(5,i,j,k,iElem)*U(5,i,j,k,iElem) &
            + U(6,i,j,k,iElem)*U(6,i,j,k,iElem)
#endif /*PP_nVar=8*/        
      ! if x, y or z is in PML region
      IF (Elem_xGP(1,i,j,k,iElem) .GE. xyzPhysicalMinMax(1) .AND. Elem_xGP(1,i,j,k,iElem) .LE. xyzPhysicalMinMax(2) .AND. &
          Elem_xGP(2,i,j,k,iElem) .GE. xyzPhysicalMinMax(3) .AND. Elem_xGP(2,i,j,k,iElem) .LE. xyzPhysicalMinMax(4) .AND. &
          Elem_xGP(3,i,j,k,iElem) .GE. xyzPhysicalMinMax(5) .AND. Elem_xGP(3,i,j,k,iElem) .LE. xyzPhysicalMinMax(6)) THEN        
          WEl_tmp  = WEl_tmp  + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * E_abs 
#if (PP_nVar==8)
          WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
#endif /*PP_nVar=8*/        
      END IF
    END DO; END DO; END DO
    WEl = WEl + WEl_tmp
#if (PP_nVar==8)
    WMag = WMag + WMag_tmp
#endif /*PP_nVar=8*/        
  END DO
ELSE
#endif /*PP_HDG*/
  DO iElem=1,nElems
    !--- Calculate and save volume of element iElem
    WEl_tmp=0. 
    WMag_tmp=0. 
    J_N(1,0:PP_N,0:PP_N,0:PP_N)=1./sJ(:,:,:,iElem)
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  ! in electromagnetische felder by henke 2011 - springer
  ! WMag = 1/(2mu) * int_V B^2 dV 

#ifdef PP_HDG
#if PP_nVar==1
      E_abs = E(1,i,j,k,iElem)*E(1,i,j,k,iElem) + E(2,i,j,k,iElem)*E(2,i,j,k,iElem) + E(3,i,j,k,iElem)*E(3,i,j,k,iElem)
#elif PP_nVar==3
      B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#else /*PP_nVar==4*/
      E_abs = E(1,i,j,k,iElem)*E(1,i,j,k,iElem) + E(2,i,j,k,iElem)*E(2,i,j,k,iElem) + E(3,i,j,k,iElem)*E(3,i,j,k,iElem)
      B_abs = B(1,i,j,k,iElem)*B(1,i,j,k,iElem) + B(2,i,j,k,iElem)*B(2,i,j,k,iElem) + B(3,i,j,k,iElem)*B(3,i,j,k,iElem)
#endif /*PP_nVar==1*/
#else
      E_abs = U(1,i,j,k,iElem)*U(1,i,j,k,iElem) + U(2,i,j,k,iElem)*U(2,i,j,k,iElem) + U(3,i,j,k,iElem)*U(3,i,j,k,iElem)
#endif /*PP_HDG*/

#if (PP_nVar==8)
      B_abs = U(4,i,j,k,iElem)*U(4,i,j,k,iElem) + U(5,i,j,k,iElem)*U(5,i,j,k,iElem) + U(6,i,j,k,iElem)*U(6,i,j,k,iElem)
#endif /*PP_nVar=8*/        

      ! if x, y or z is in PML region
#ifdef PP_HDG
#if PP_nVar==3
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
#elif PP_nVar==4
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
#endif /*PP_nVar==3*/
#endif /*PP_HDG*/
      WEl_tmp  = WEl_tmp  + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * E_abs 
#if (PP_nVar==8)
      WMag_tmp = WMag_tmp + wGP(i)*wGP(j)*wGP(k) * J_N(1,i,j,k) * B_abs
#endif /*PP_nVar=8*/        
    END DO; END DO; END DO
    WEl = WEl + WEl_tmp
#if (PP_nVar==8)
    WMag = WMag + WMag_tmp
#endif /*PP_nVar=8*/        
  END DO
#ifndef PP_HDG
END IF ! noPML
#endif /*PP_HDG*/

WEl = WEl * eps0 * 0.5 
WMag = WMag * smu0 * 0.5

#ifdef MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,WEl  , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, IERROR)
  CALL MPI_REDUCE(MPI_IN_PLACE,WMag , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, IERROR)
ELSE
  CALL MPI_REDUCE(WEl         ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, IERROR)
  CALL MPI_REDUCE(WMag        ,RD   , 1 , MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, IERROR)
END IF
#endif /*MPI*/

END SUBROUTINE CalcPotentialEnergy


END MODULE MOD_AnalyzeField
