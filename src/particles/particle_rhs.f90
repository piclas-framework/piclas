#include "boltzplatz.h"

MODULE MOD_part_RHS                                                                                
!===================================================================================================================================
! 
!===================================================================================================================================
IMPLICIT NONE                                                                                    
PRIVATE                                                                                          
!----------------------------------------------------------------------------------------------------------------------------------
PUBLIC            :: CalcPartRHS                                                                 
!----------------------------------------------------------------------------------------------------------------------------------
                                                                                                   
CONTAINS                                                                                           
                                                                                                   
SUBROUTINE CalcPartRHS()                                                                           
!===================================================================================================================================
!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Particle_Vars,          ONLY : PDM, PartState, Pt, Species, PartSpecies, PartLorentzType
USE MOD_Equation_Vars,          ONLY : c2_inv, c ,c2
USE MOD_TimeDisc_Vars,          ONLY : dt
USE MOD_PICInterpolation_Vars,  ONLY : FieldAtParticle                                           
USE MOD_part_MPI_Vars,          ONLY :  PMPIVAR
!----------------------------------------------------------------------------------------------------------------------------------
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE                                                                                    
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLE 
INTEGER                          :: counter1,i,n                                                 
REAL                             :: E(1:3),B(1:3)                                                
REAL                             :: qmt, gamma, velosq, v_norm_x, v_norm_y, v_norm_z,v_abs, Pt_sq
REAL                             :: ax, ay, az, bx, by, bz, dx, dy, snx, sny, snz
REAL                             :: gamma1,gamma2,gamma3, v1s,v2s,v3s, Vinv(3,3), v1,v2,v3, normfac
!===================================================================================================================================

! Lorentzforce
Pt(1:PDM%ParticleVecLength,:)=0.
SELECT CASE(PartLorentzType)
  CASE(1) ! default
    ! constant Lorentz factor over time step
    DO i = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(i)) THEN
      ! Calculation of relativistic Factor: m_rel = m0 * 1/sqrt(1-|v^2/c^2|)

        velosq = PartState(i,4) * PartState(i,4) &
                 + PartState(i,5) * PartState(i,5) &
                 + PartState(i,6) * PartState(i,6)  
      ! MPF in ChargeIC and MassIC cancels out.
        ! old method
        gamma = (SQRT(1.0 - velosq * c2_inv))
        qmt = Species(PartSpecies(i))%ChargeIC/Species(PartSpecies(i))%MassIC * gamma
        E(1:3) = FieldAtParticle(i,1:3) * qmt
#if (PP_nVar==8)
        B(1:3) = FieldAtParticle(i,4:6) * qmt
#endif
      ! Calc Lorentz forces in x, y, z direction:
#if (PP_nVar==8)
        Pt(i,1) = E(1) + PartState(i,5) * B(3) - PartState(i,6) * B(2)
        Pt(i,2) = E(2) + PartState(i,6) * B(1) - PartState(i,4) * B(3)
        Pt(i,3) = E(3) + PartState(i,4) * B(2) - PartState(i,5) * B(1)
#else
        Pt(i,1) = E(1) 
        Pt(i,2) = E(2) 
        Pt(i,3) = E(3) 
#endif
      END IF
    END DO
  CASE(2)
  ! exact Lorentz-Pusher
  ! prevent particles from acceleration above speed of light
    DO i = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(i)) THEN
    ! Calculation of relativistic Factor: m_rel = m0 * 1/sqrt(1-|v^2/c^2|)
        velosq = PartState(i,4) * PartState(i,4) &
               + PartState(i,5) * PartState(i,5) &
               + PartState(i,6) * PartState(i,6)  
          ! MPF in ChargeIC and MassIC cancels out.
        qmt = Species(PartSpecies(i))%ChargeIC/Species(PartSpecies(i))%MassIC
        E(1:3) = FieldAtParticle(i,1:3) * qmt
#if (PP_nVar==8)
        B(1:3) = FieldAtParticle(i,4:6) * qmt
#endif
      ! Calc Lorentz forces in x, y, z direction:
#if (PP_nVar==8)
        Pt(i,1) = E(1) + PartState(i,5) * B(3) - PartState(i,6) * B(2)
        Pt(i,2) = E(2) + PartState(i,6) * B(1) - PartState(i,4) * B(3)
        Pt(i,3) = E(3) + PartState(i,4) * B(2) - PartState(i,5) * B(1)
#else
        Pt(i,1) = E(1) 
        Pt(i,2) = E(2) 
        Pt(i,3) = E(3) 
#endif

        gamma = 1/sqrt(1.0 - velosq * c2_inv)
        bx = Pt(i,1) *dt + gamma * PartState(i,4)
        snx = sign(1.0,bx)
        bx = bx*bx*c2_inv
        bx = bx/(1+bx)

        by = Pt(i,2) *dt + gamma * PartState(i,5)
        sny = sign(1.0,by)
        by = by*by*c2_inv
        by = by/(1+by)

        bz = Pt(i,3) *dt + gamma * PartState(i,6)
        snz = sign(1.0,bz)
        bz = bz*bz*c2_inv
        bz = bz/(1+bz)

        dx = (bx-bx*bz)/(1-bx*bz)
        dy = (by-by*bz)/(1-by*bz)
          
        ax = (dx-dx*dy)/(1-dx*dy)
        ay = (dy-dy*ax)
        az = (bz*(1-ax-ay))

        Pt(i,1) = (snx * sqrt(ax)*c - PartState(i,4)) / dt
        Pt(i,2) = (sny * sqrt(ay)*c - PartState(i,5)) / dt
        Pt(i,3) = (snz * sqrt(az)*c - PartState(i,6)) / dt

!       Pt_sq = Pt(i,1)*Pt(i,1)+Pt(i,2)*Pt(i,2)+Pt(i,3)*Pt(i,3)
!       IF(Pt_sq.EQ.0.0) THEN
!         v_norm_x = 0
!         v_norm_y = 0
!         v_norm_z = 0
!       ELSE
!         v_norm_x = Pt(i,1) / SQRT(Pt_sq)
!         v_norm_y = Pt(i,2) / SQRT(Pt_sq)
!         v_norm_z = Pt(i,3) / SQRT(Pt_sq)
!       END IF
!
!        v_abs = (SQRT(Pt_sq) * dt + SQRT(velosq)/SQRT(1.0-velosq * c2_inv)) &
!                / SQRT(1.0 + (SQRT(Pt_sq) * dt + SQRT(velosq)/SQRT(1.0-velosq * c2_inv))**2 * c2_inv)
!
!       if (IsNAN(v_abs))print*, "Geschw: ", v_abs, i, PMPIVAR%iProc, SQRT(velosq), ASIN(SQRT(velosq)/c), c, SQRT(velosq)/c
!
!        Pt(i,1) = (v_abs - SQRT(velosq)) / dt * v_norm_x
!        Pt(i,2) = (v_abs - SQRT(velosq)) / dt * v_norm_y
!        Pt(i,3) = (v_abs - SQRT(velosq)) / dt * v_norm_z
      END IF
    END DO

  CASE(3)
    ! derivation of relativistic equation of motion
    DO i = 1,PDM%ParticleVecLength
      IF (PDM%ParticleInside(i)) THEN
        ! Calculation of relativistic Factor: m_rel = m0 * 1/sqrt(1-|v^2/c^2|)

        ! required helps
        v1  = PartState(i,4)
        v2  = PartState(i,5)
        v3  = PartState(i,6)

        v1s = PartState(i,4) * PartState(i,4)
        v2s = PartState(i,5) * PartState(i,5)
        v3s = PartState(i,6) * PartState(i,6)
        velosq = v1s+v2s+v3s

        gamma=SQRT(1.0 - velosq*c2_inv)
        gamma2=gamma1*gamma1
        gamma3=gamma2*gamma1
        normfac=1.0/(c2*gamma2 + velosq)
        ! define inverted matrix
        Vinv(1,1) = (c2*gamma3 + (v2s+v3s)*gamma)*normfac
        Vinv(1,2) =-(v1*v2*gamma)*normfac
        Vinv(1,3) =-(v1*v3*gamma)*normfac
        Vinv(2,1) = Vinv(1,2)
        Vinv(2,2) = (c2*gamma3 + (v1s+v3s)*gamma)*normfac
        Vinv(2,3) =-(v2*v3*gamma)*normfac
        Vinv(3,1) = Vinv(1,3)
        Vinv(3,2) = Vinv(2,3)
        Vinv(3,3) = (c2*gamma3 + (v1s+v2s)*gamma)*normfac

!        gamma1= 1./SQRT(1.0 - velosq * c2_inv)
!        gamma2=gamma1*gamma1
!        gamma3=gamma2*gamma1
!        ! define inverted matrix
!        Vinv(1,1) = ((v3s+v2s)*gamma2+c2)/(velosq*gamma3+c2*gamma1)
!        Vinv(1,2) =-v1*v2*gamma1/(velosq*gamma2+c2)
!        Vinv(1,2) =-v1*v3*gamma1/(velosq*gamma2+c2)
!        Vinv(2,1) = Vinv(1,2)
!        Vinv(2,2) = ((v3s+v1s)*gamma2+c2)/(velosq*gamma3+c2*gamma1)
!        Vinv(2,3) =-v2*v3*gamma1/(velosq*gamma2+c2)
!        Vinv(3,1) = Vinv(1,3)
!        Vinv(3,2) = Vinv(2,3)
!        Vinv(3,3) =((v1s+v2s)*gamma2+c2)/(velosq*gamma3+c2*gamma1)

        qmt = Species(PartSpecies(i))%ChargeIC/Species(PartSpecies(i))%MassIC

        E(1:3) = FieldAtParticle(i,1:3) * qmt
#if (PP_nVar==8)
        B(1:3) = FieldAtParticle(i,4:6) * qmt
#endif
      ! Calc Lorentz forces in x, y, z direction:
#if (PP_nVar==8)
        Pt(i,1) = E(1) + PartState(i,5) * B(3) - PartState(i,6) * B(2)
        Pt(i,2) = E(2) + PartState(i,6) * B(1) - PartState(i,4) * B(3)
        Pt(i,3) = E(3) + PartState(i,4) * B(2) - PartState(i,5) * B(1)
#else
        Pt(i,1) = E(1) 
        Pt(i,2) = E(2) 
        Pt(i,3) = E(3) 
#endif
        Pt(i,1:3) = MATMUL(Vinv,Pt(i,1:3))

      END IF
    END DO
  CASE DEFAULT
    CALL abort(__STAMP__,&
      'This Type of Lorentz-force calculation is not implemented:.',PartLorentzType,999.)
END SELECT

END SUBROUTINE CalcPartRHS

END MODULE MOD_part_RHS
