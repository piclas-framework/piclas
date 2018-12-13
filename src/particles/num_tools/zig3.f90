!==================================================================================================================================
!
! File from: https://jblevins.org/mirror/amiller/
! -----------------------------------------------
! ziggurat.f90 George Marsaglia's functions for generating random samples from the uniform, normal and exponential distributions.
! Translated from C.2
!
! 2018: Modified for user-defined seed !!!
!==================================================================================================================================
! Marsaglia & Tsang generator for random normals & random exponentials.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001
MODULE Ziggurat
   IMPLICIT NONE

   PRIVATE

   INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )
   REAL(DP), PARAMETER  ::  m1=2147483648.0_DP,   m2=2147483648.0_DP,      &
                            half=0.5_DP
   REAL(DP)             ::  dn=3.442619855899_DP, tn=3.442619855899_DP,    &
                            vn=0.00991256303526217_DP,                     &
                            q,                    de=7.697117470131487_DP, &
                            te=7.697117470131487_DP,                       &
                            ve=0.003949659822581572_DP
   INTEGER,  SAVE       ::  iz, jz, jsr, kn(0:127),              &
                            ke(0:255), hz
   REAL(DP), SAVE       ::  wn(0:127), fn(0:127), we(0:255), fe(0:255)
   LOGICAL,  SAVE       ::  initialized=.FALSE.

   PUBLIC  :: zigset, shr3, uni, rnor, rexp


CONTAINS


SUBROUTINE zigset( ) !jsrseed )
   USE MOD_Particle_Vars,      ONLY : seeds

   !INTEGER, INTENT(IN)  :: jsrseed

   INTEGER  :: i

   !  Set the seed
   jsr = seeds(1) !jsrseed

   !  Tables for RNOR
   q = vn*EXP(half*dn*dn)
   kn(0) = INT((dn/q)*m1)
   kn(1) = 0
   wn(0) = q/m1
   wn(127) = dn/m1
   fn(0) = 1.0_DP
   fn(127) = EXP( -half*dn*dn )
   DO  i = 126, 1, -1
      dn = SQRT( -2.0_DP * LOG( vn/dn + EXP( -half*dn*dn ) ) )
      kn(i+1) = INT((dn/tn)*m1)
      tn = dn
      fn(i) = EXP(-half*dn*dn)
      wn(i) = dn/m1
   END DO

   !  Tables for REXP
   q = ve*EXP( de )
   ke(0) = INT((de/q)*m2)
   ke(1) = 0
   we(0) = q/m2
   we(255) = de/m2
   fe(0) = 1.0_DP
   fe(255) = EXP( -de )
   DO  i = 254, 1, -1
      de = -LOG( ve/de + EXP( -de ) )
      ke(i+1) = INT(m2 * (de/te))
      te = de
      fe(i) = EXP( -de )
      we(i) = de/m2
   END DO
   initialized = .TRUE.
   RETURN
END SUBROUTINE zigset



!  Generate random 32-bit integers
FUNCTION shr3( ) RESULT( ival )
   INTEGER  ::  ival

   jz = jsr
   jsr = IEOR( jsr, ISHFT( jsr,  13 ) )
   jsr = IEOR( jsr, ISHFT( jsr, -17 ) )
   jsr = IEOR( jsr, ISHFT( jsr,   5 ) )
   ival = jz + jsr
   RETURN
END FUNCTION shr3



!  Generate uniformly distributed random numbers
FUNCTION uni( ) RESULT( fn_val )
   REAL(DP)  ::  fn_val

   fn_val = half + 0.2328306e-9_DP * shr3( )
   RETURN
END FUNCTION uni



!  Generate random normals
FUNCTION rnor( ) RESULT( fn_val )
   REAL(DP)             ::  fn_val

   REAL(DP), PARAMETER  ::  r = 3.442620_DP
   REAL(DP)             ::  x, y

   IF( .NOT. initialized ) CALL zigset( ) !jsr )
   hz = shr3( )
   iz = IAND( hz, 127 )
   IF( ABS( hz ) < kn(iz) ) THEN
      fn_val = hz * wn(iz)
   ELSE
      DO
         IF( iz == 0 ) THEN
            DO
               x = -0.2904764_DP* LOG( uni( ) )
               y = -LOG( uni( ) )
               IF( y+y >= x*x ) EXIT
            END DO
            fn_val = r+x
            IF( hz <= 0 ) fn_val = -fn_val
            RETURN
         END IF
         x = hz * wn(iz)
         IF( fn(iz) + uni( )*(fn(iz-1)-fn(iz)) < EXP(-half*x*x) ) THEN
            fn_val = x
            RETURN
         END IF
         hz = shr3( )
         iz = IAND( hz, 127 )
         IF( ABS( hz ) < kn(iz) ) THEN
            fn_val = hz * wn(iz)
            RETURN
         END IF
      END DO
   END IF
   RETURN
END FUNCTION rnor



!  Generate random exponentials
FUNCTION rexp( ) RESULT( fn_val )
   REAL(DP)  ::  fn_val

   REAL(DP)  ::  x

   IF( .NOT. initialized ) CALL Zigset( ) !jsr )
   jz = shr3( )
   iz = IAND( jz, 255 )
   IF( ABS( jz ) < ke(iz) ) THEN
      fn_val = ABS(jz) * we(iz)
      RETURN
   END IF
   DO
      IF( iz == 0 ) THEN
         fn_val = 7.69711 - LOG( uni( ) )
         RETURN
      END IF
      x = ABS( jz ) * we(iz)
      IF( fe(iz) + uni( )*(fe(iz-1) - fe(iz)) < EXP( -x ) ) THEN
         fn_val = x
         RETURN
      END IF
      jz = shr3( )
      iz = IAND( jz, 255 )
      IF( ABS( jz ) < ke(iz) ) THEN
         fn_val = ABS( jz ) * we(iz)
         RETURN
      END IF
   END DO
   RETURN
END FUNCTION rexp

END MODULE ziggurat
