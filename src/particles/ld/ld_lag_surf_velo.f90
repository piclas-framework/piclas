MODULE MOD_LD_lag_velo

!===================================================================================================================================
! module for determination of Lagrangian cell velocity
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

PUBLIC :: CalcSurfLagVelo
!===================================================================================================================================

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE CalcSurfLagVelo

USE MOD_LD_Vars
USE MOD_Mesh_Vars,             ONLY : nElems, SideToElem, BC, ElemToSide
USE MOD_Particle_Vars,         ONLY : PartBound
USE MOD_TimeDisc_Vars,         ONLY : iter
#ifdef MPI
USE MOD_Mesh_Vars,             ONLY : nBCSides, nInnerSides
#endif

!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
! Local variable declaration                                                                       !
  INTEGER           :: iElem, trinum, Elem2, IterForSecant, iLocSide, SideID, iLocSide2
  REAL              :: Velo1(3), Velo2(3)
  REAL              :: Beta1, Beta2, Dens1, Dens2
  REAL              :: VeloDiff1_old, VeloDiff2_old, VeloDiff1_new, VeloDiff2_new
  REAL              :: NVec(3)
  REAL              :: kon1, kon2, vLAG_old, vLAG_new, vLAG, VeloDir1, VeloDir2
  REAL              :: G_old, G_new
  REAL, PARAMETER   :: PI=3.14159265358979323846_8
  LOGICAL           :: IsStationary
!--------------------------------------------------------------------------------------------------!
  DO iElem = 1, nElems
    Velo1 = BulkValues(iElem)%CellV
    Beta1 = BulkValues(iElem)%Beta
    Dens1 = BulkValues(iElem)%MassDens
    DO iLocSide = 1, 6
      Elem2 = 0
      SideID = ElemToSide(1,iLocSide,iElem)
      IF (.NOT.IsDoneLagVelo(SideID)) THEN
        IsDoneLagVelo(SideID) = .TRUE.
        IsStationary = .FALSE.
#ifdef MPI
        IF (SideID.GT.nBCSides+nInnerSides) THEN ! it must be a MPI Side
          Elem2 = -1
          Velo2(1) = MPINeighborBulkVal(SideID,1)
          Velo2(2) = MPINeighborBulkVal(SideID,2) 
          Velo2(3) = MPINeighborBulkVal(SideID,3)
          Beta2 = MPINeighborBulkVal(SideID,4)
          Dens2 = MPINeighborBulkVal(SideID,5)
        ELSE
#endif
        IF (SideToElem(1,SideID) .EQ. iElem) THEN
          IF (SideToElem(2,SideID) .LT. 1) THEN ! it must be a BC
            IF (PartBound%Map(BC(SideID)).EQ.PartBound%OpenBC) THEN ! open => copy cell values
              Elem2 = -1
              IF (iter.GE. 1) THEN
                IF (PartBound%AmbientCondition(BC(SideID))) THEN  
                  Velo2(1) = PartBound%AmbientVelo(1,BC(SideID))
                  Velo2(2) = PartBound%AmbientVelo(2,BC(SideID))
                  Velo2(3) = PartBound%AmbientVelo(3,BC(SideID))
                  Beta2 = PartBound%AmbientBeta(BC(SideID))
                  Dens2 = PartBound%AmbientDens(BC(SideID))
                ELSE  
                  BulkValuesOpenBC(iElem)%CellV = (1 - LD_RelaxationFak) * BulkValuesOpenBC(iElem)%CellV &
                                                + LD_RelaxationFak * BulkValues(iElem)%CellV
                  BulkValuesOpenBC(iElem)%Beta = (1 - LD_RelaxationFak) * BulkValuesOpenBC(iElem)%Beta &
                                               + LD_RelaxationFak * BulkValues(iElem)%Beta
                  BulkValuesOpenBC(iElem)%MassDens = (1 - LD_RelaxationFak) * BulkValuesOpenBC(iElem)%MassDens & 
                                                   + LD_RelaxationFak * BulkValues(iElem)%MassDens
                  Velo2 = BulkValuesOpenBC(iElem)%CellV
                  Beta2 = BulkValuesOpenBC(iElem)%Beta
                  Dens2 = BulkValuesOpenBC(iElem)%MassDens
                END IF
              ELSE
                IF (PartBound%AmbientCondition(BC(SideID))) THEN  
                  Velo2(1) = PartBound%AmbientVelo(1,BC(SideID))
                  Velo2(2) = PartBound%AmbientVelo(2,BC(SideID))
                  Velo2(3) = PartBound%AmbientVelo(3,BC(SideID))
                  Beta2 = PartBound%AmbientBeta(BC(SideID))
                  Dens2 = PartBound%AmbientDens(BC(SideID))
                ELSE  
                  BulkValuesOpenBC(iElem)%CellV    = BulkValues(iElem)%CellV
                  BulkValuesOpenBC(iElem)%Beta     = BulkValues(iElem)%Beta
                  BulkValuesOpenBC(iElem)%MassDens = BulkValues(iElem)%MassDens
                  Velo2 = BulkValuesOpenBC(iElem)%CellV
                  Beta2 = BulkValuesOpenBC(iElem)%Beta
                  Dens2 = BulkValuesOpenBC(iElem)%MassDens
                END IF
              END IF
            ELSE IF (PartBound%Map(BC(SideID)).EQ.PartBound%ReflectiveBC) THEN
              IsStationary = .TRUE.
            ELSE
              WRITE(*,*)'=============================================='
              WRITE(*,*)' ERROR in PartBound%Map(BC(SideID))'
              STOP
            END IF
          ELSE
            IF (SideToElem(1,SideID).EQ.SideToElem(2,SideID)) THEN ! one periodic cell
              Elem2 = iElem
              SELECT CASE(iLocSide)
                CASE (1)
                  iLocSide2 = 6
                CASE (2)
                  iLocSide2 = 4
                CASE (3)
                  iLocSide2 = 5
                CASE DEFAULT
                  WRITE(*,*)'=============================================='
                  WRITE(*,*)' ERROR in LocSides for periodic Element'
                  STOP
              END SELECT
              Velo2 = BulkValues(iElem)%CellV
              Beta2 = BulkValues(iElem)%Beta
              Dens2 = BulkValues(iElem)%MassDens            
            ELSE
              Elem2 = SideToElem(2,SideID)
              iLocSide2 = SideToElem(4,SideID)
              Velo2 = BulkValues(Elem2)%CellV
              Beta2 = BulkValues(Elem2)%Beta
              Dens2 = BulkValues(Elem2)%MassDens
            END IF
          END IF
        ELSE         
          IF (SideToElem(1,SideID) .LT. 1) THEN ! it must be a BC
            IF (PartBound%Map(BC(SideID)).EQ.PartBound%OpenBC) THEN ! open => copy cell values
              Elem2 = -1
              IF (PartBound%AmbientCondition(BC(SideID))) THEN  
                Velo2(1) = PartBound%AmbientVelo(1,BC(SideID))
                Velo2(2) = PartBound%AmbientVelo(2,BC(SideID))
                Velo2(3) = PartBound%AmbientVelo(3,BC(SideID))
                Beta2 = PartBound%AmbientBeta(BC(SideID))
                Dens2 = PartBound%AmbientDens(BC(SideID))
              ELSE 
                Velo2 = BulkValues(iElem)%CellV
                Beta2 = BulkValues(iElem)%Beta
                Dens2 = BulkValues(iElem)%MassDens
              END IF
            ELSE IF (PartBound%Map(BC(SideID)).EQ.PartBound%ReflectiveBC) THEN
              IsStationary = .TRUE.
            ELSE
              WRITE(*,*)'=============================================='
              WRITE(*,*)' ERROR in PartBound%Map(BC(SideID))'
              STOP
            END IF
          ELSE
            Elem2 = SideToElem(1,SideID)
            iLocSide2 = SideToElem(3,SideID)
            Velo2 = BulkValues(Elem2)%CellV
            Beta2 = BulkValues(Elem2)%Beta
            Dens2 = BulkValues(Elem2)%MassDens
          END IF
        END IF
#ifdef MPI
        END IF
#endif
        IF (.NOT.IsStationary) THEN
          DO trinum=1, 3  ! third loop for mean surf value
            IF (trinum.EQ.3) THEN
              NVec = MeanSurfValues(iLocSide, iElem)%MeanNormVec
            ELSE
              NVec = SurfLagValues(iLocSide, iElem, trinum)%LagNormVec
            END IF
            kon1 = Dens1 / Beta1**2
            kon2 = Dens2 / Beta2**2
            VeloDir1 = Velo1(1) * NVec(1) &
                     + Velo1(2) * NVec(2) &  
                     + Velo1(3) * NVec(3) 
            VeloDir2 = Velo2(1) * NVec(1) &
                     + Velo2(2) * NVec(2) &
                     + Velo2(3) * NVec(3) 
            IF (trinum.EQ.3) THEN
              vLAG_old = MeanSurfValues(iLocSide, iElem)%MeanLagVelo   ! 1st guess, former velocity, [m/s]
              vLAG_new = MeanSurfValues(iLocSide, iElem)%MeanLagVelo &
                       + LD_SecantMeth%Guess                     ! 2nd guess, plus user defined value, [m/s], (default 10 m/s)
            ELSE
              vLAG_old = SurfLagValues(iLocSide, iElem, trinum)%LagVelo   ! 1st guess, former velocity, [m/s]
              vLAG_new = SurfLagValues(iLocSide, iElem, trinum)%LagVelo &
                       + LD_SecantMeth%Guess                     ! 2nd guess, plus user defined value, [m/s], (default 10 m/s)
            END IF
            IterForSecant = 0                                           ! reset iteration counter
      !
      !.... find root of G-function ==> Lag. Surface Velocity
      !
            DO WHILE (ABS(vLAG_old - vLAG_new) .GT. LD_SecantMeth%Accuracy)  ! iteration loop
              IF (IterForSecant .GT. LD_SecantMeth%MaxIter) THEN
                WRITE(*,*)'=============================================='
                WRITE(*,*)' Max. number of iterations for LAGRANGian cell'
                WRITE(*,*)' exceeded. Change something and restart run...'              
                PRINT*,iLocSide, iElem
                PRINT*,Velo1, Velo2
                PRINT*,Beta1, Beta2
                PRINT*,Dens1, Dens2
                STOP
              END IF
              IterForSecant = IterForSecant + 1                                   ! increase local iteration counter
              VeloDiff1_new =  Beta1 * ( VeloDir1 - vLAG_new )
              VeloDiff2_new = -Beta2 * ( VeloDir2 - vLAG_new ) !
              VeloDiff1_old =  Beta1 * ( VeloDir1 - vLAG_old ) 
              VeloDiff2_old = -Beta2 * ( VeloDir2 - vLAG_old ) !
              G_old  = kon1 * ( VeloDiff1_old * EXP(-VeloDiff1_old**2) + SQRT(PI) &
                     * ( 1. + ERF(VeloDiff1_old) ) * ( 0.5 + VeloDiff1_old**2 ) ) &
                     - kon2 * ( VeloDiff2_old * EXP(-VeloDiff2_old**2) + SQRT(PI) &
                     * ( 1. + ERF(VeloDiff2_old) ) * ( 0.5 + VeloDiff2_old**2) )
              G_new  = kon1 * ( VeloDiff1_new * EXP(-VeloDiff1_new**2) + SQRT(PI) &
                     * ( 1. + ERF(VeloDiff1_new) ) * ( 0.5 + VeloDiff1_new**2 ) ) &
                     - kon2 * ( VeloDiff2_new * EXP(-VeloDiff2_new**2) + SQRT(PI) &
                     * ( 1. + ERF(VeloDiff2_new) ) * ( 0.5 + VeloDiff2_new**2) ) 
              IF ( ( G_new - G_old ) .NE. 0.0 ) THEN
                vLAG  = vLAG_new + ( vLAG_old - vLAG_new ) * ( G_new / ( G_new - G_old ) )
              ELSE
                vLAG  = vLAG_new
              END IF
              vLAG_old = vLAG_new                                        ! update start values
              vLAG_new = vLAG
            END DO                                                      ! end of iteration loop
            IF (trinum.EQ.3) THEN
                MeanSurfValues(iLocSide, iElem)%MeanLagVelo = vLAG
                IF (Elem2.GT. 0) THEN
                  MeanSurfValues(iLocSide2, Elem2)%MeanLagVelo =  (-1.) * vLAG
                END IF
            ELSE
                SurfLagValues(iLocSide, iElem,trinum)%LagVelo = vLAG
                IF (Elem2.GT. 0) SurfLagValues(iLocSide2, Elem2,trinum)%LagVelo =  (-1.) * vLAG
            END IF
          END DO  ! end of trinum
        ELSE ! Side is stationary
          MeanSurfValues(iLocSide, iElem)%MeanLagVelo = 0.0
          SurfLagValues(iLocSide, iElem,1)%LagVelo = 0.0
          SurfLagValues(iLocSide, iElem,2)%LagVelo = 0.0
        END IF
      END IF ! end if done
    END DO ! end loop over ilocsides
  END DO ! end loop over elements

END SUBROUTINE CalcSurfLagVelo


END MODULE MOD_LD_lag_velo
