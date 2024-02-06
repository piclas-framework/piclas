!=================================================================================================================================
! Copyright (c) 2010-2021  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "piclas.h"

!==================================================================================================================================
!> Routines providing distribution function management for the exponential integration method with discrete velocities
!==================================================================================================================================
MODULE MOD_IncBeta
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC:: betain
!==================================================================================================================================

CONTAINS

function alogam ( x, ifault )

    !*****************************************************************************80
    !
    !! ALOGAM computes the logarithm of the Gamma function.
    !
    !  Modified:
    !
    !    28 March 1999
    !
    !  Author:
    !
    !    Malcolm Pike, David Hill.
    !    FORTRAN90 version by John Burkardt
    !
    !  Reference:
    !
    !    Malcolm Pike, David Hill,
    !    Algorithm 291:
    !    Logarithm of Gamma Function,
    !    Communications of the ACM,
    !    Volume 9, Number 9, September 1966, page 684.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the Gamma function.
    !    X should be greater than 0.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error flag.
    !    0, no error.
    !    1, X <= 0.
    !
    !    Output, real ( kind = 8 ) ALOGAM, the logarithm of the Gamma
    !    function of X.
    !
      implicit none

      real ( kind = 8 ) alogam
      real ( kind = 8 ) f
      integer ( kind = 4 ) ifault
      real ( kind = 8 ) x
      real ( kind = 8 ) y
      real ( kind = 8 ) z

      if ( x <= 0.0D+00 ) then
        ifault = 1
        alogam = 0.0D+00
        return
      end if

      ifault = 0
      y = x

      if ( x < 7.0D+00 ) then

        f = 1.0D+00
        z = y

        do while ( z < 7.0D+00 )
          f = f * z
          z = z + 1.0D+00
        end do

        y = z
        f = - log ( f )

      else

        f = 0.0D+00

      end if

      z = 1.0D+00 / y / y

      alogam = f + ( y - 0.5D+00 ) * log ( y ) - y &
        + 0.918938533204673D+00 + &
        ((( &
        - 0.000595238095238D+00   * z &
        + 0.000793650793651D+00 ) * z &
        - 0.002777777777778D+00 ) * z &
        + 0.083333333333333D+00 ) / y

      return
    end
    subroutine beta_inc_values ( n_data, a, b, x, fx )

    !*****************************************************************************80
    !
    !! BETA_INC_VALUES returns some values of the incomplete Beta function.
    !
    !  Discussion:
    !
    !    The incomplete Beta function may be written
    !
    !      BETA_INC(A,B,X) = integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
    !                      / integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
    !
    !    Thus,
    !
    !      BETA_INC(A,B,0.0) = 0.0
    !      BETA_INC(A,B,1.0) = 1.0
    !
    !    The incomplete Beta function is also sometimes called the
    !    "modified" Beta function, or the "normalized" Beta function
    !    or the Beta CDF (cumulative density function).
    !
    !    In Mathematica, the function can be evaluated by:
    !
    !      BETA[X,A,B] / BETA[A,B]
    !
    !    The function can also be evaluated by using the Statistics package:
    !
    !      Needs["Statistics`ContinuousDistributions`"]
    !      dist = BetaDistribution [ a, b ]
    !      CDF [ dist, x ]
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 April 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    National Bureau of Standards, 1964,
    !    ISBN: 0-486-61272-4,
    !    LC: QA47.A34.
    !
    !    Karl Pearson,
    !    Tables of the Incomplete Beta Function,
    !    Cambridge University Press, 1968,
    !    LC: QA351.P38.
    !
    !    Stephen Wolfram,
    !    The Mathematica Book,
    !    Fourth Edition,
    !    Cambridge University Press, 1999,
    !    ISBN: 0-521-64314-7,
    !    LC: QA76.95.W65.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
    !    before the first call.  On each call, the routine increments N_DATA by 1,
    !    and returns the corresponding data; when there is no more data, the
    !    output value of N_DATA will be 0 again.
    !
    !    Output, real ( kind = 8 ) A, B, the parameters of the function.
    !
    !    Output, real ( kind = 8 ) X, the argument of the function.
    !
    !    Output, real ( kind = 8 ) FX, the value of the function.
    !
      implicit none

      integer ( kind = 4 ), parameter :: n_max = 45

      real ( kind = 8 ) a
      real ( kind = 8 ), save, dimension ( n_max ) :: a_vec = (/ &
         0.5D+00, &
         0.5D+00, &
         0.5D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         5.5D+00, &
        10.0D+00, &
        10.0D+00, &
        10.0D+00, &
        10.0D+00, &
        20.0D+00, &
        20.0D+00, &
        20.0D+00, &
        20.0D+00, &
        20.0D+00, &
        30.0D+00, &
        30.0D+00, &
        40.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         1.0D+00, &
         2.0D+00, &
         3.0D+00, &
         4.0D+00, &
         5.0D+00, &
         1.30625D+00, &
         1.30625D+00, &
         1.30625D+00 /)
      real ( kind = 8 ) b
      real ( kind = 8 ), save, dimension ( n_max ) :: b_vec = (/ &
         0.5D+00, &
         0.5D+00, &
         0.5D+00, &
         0.5D+00, &
         0.5D+00, &
         0.5D+00, &
         0.5D+00, &
         1.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         5.0D+00, &
         0.5D+00, &
         5.0D+00, &
         5.0D+00, &
        10.0D+00, &
         5.0D+00, &
        10.0D+00, &
        10.0D+00, &
        20.0D+00, &
        20.0D+00, &
        10.0D+00, &
        10.0D+00, &
        20.0D+00, &
         0.5D+00, &
         0.5D+00, &
         0.5D+00, &
         0.5D+00, &
         2.0D+00, &
         3.0D+00, &
         4.0D+00, &
         5.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
         2.0D+00, &
        11.7562D+00, &
        11.7562D+00, &
        11.7562D+00 /)
      real ( kind = 8 ) fx
      real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.6376856085851985D-01, &
        0.2048327646991335D+00, &
        0.1000000000000000D+01, &
        0.0000000000000000D+00, &
        0.5012562893380045D-02, &
        0.5131670194948620D-01, &
        0.2928932188134525D+00, &
        0.5000000000000000D+00, &
        0.2800000000000000D-01, &
        0.1040000000000000D+00, &
        0.2160000000000000D+00, &
        0.3520000000000000D+00, &
        0.5000000000000000D+00, &
        0.6480000000000000D+00, &
        0.7840000000000000D+00, &
        0.8960000000000000D+00, &
        0.9720000000000000D+00, &
        0.4361908850559777D+00, &
        0.1516409096347099D+00, &
        0.8978271484375000D-01, &
        0.1000000000000000D+01, &
        0.5000000000000000D+00, &
        0.4598773297575791D+00, &
        0.2146816102371739D+00, &
        0.9507364826957875D+00, &
        0.5000000000000000D+00, &
        0.8979413687105918D+00, &
        0.2241297491808366D+00, &
        0.7586405487192086D+00, &
        0.7001783247477069D+00, &
        0.5131670194948620D-01, &
        0.1055728090000841D+00, &
        0.1633399734659245D+00, &
        0.2254033307585166D+00, &
        0.3600000000000000D+00, &
        0.4880000000000000D+00, &
        0.5904000000000000D+00, &
        0.6723200000000000D+00, &
        0.2160000000000000D+00, &
        0.8370000000000000D-01, &
        0.3078000000000000D-01, &
        0.1093500000000000D-01, &
        0.918884684620518D+00, &
        0.21052977489419D+00, &
        0.1824130512500673D+00 /)
      integer ( kind = 4 ) n_data
      real ( kind = 8 ) x
      real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
        0.01D+00, &
        0.10D+00, &
        1.00D+00, &
        0.00D+00, &
        0.01D+00, &
        0.10D+00, &
        0.50D+00, &
        0.50D+00, &
        0.10D+00, &
        0.20D+00, &
        0.30D+00, &
        0.40D+00, &
        0.50D+00, &
        0.60D+00, &
        0.70D+00, &
        0.80D+00, &
        0.90D+00, &
        0.50D+00, &
        0.90D+00, &
        0.50D+00, &
        1.00D+00, &
        0.50D+00, &
        0.80D+00, &
        0.60D+00, &
        0.80D+00, &
        0.50D+00, &
        0.60D+00, &
        0.70D+00, &
        0.80D+00, &
        0.70D+00, &
        0.10D+00, &
        0.20D+00, &
        0.30D+00, &
        0.40D+00, &
        0.20D+00, &
        0.20D+00, &
        0.20D+00, &
        0.20D+00, &
        0.30D+00, &
        0.30D+00, &
        0.30D+00, &
        0.30D+00, &
        0.225609D+00, &
        0.0335568D+00, &
        0.0295222D+00 /)

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
    end
    function betain ( x, p, q, beta, ifault )

    !*****************************************************************************80
    !
    !! BETAIN computes the incomplete Beta function ratio.
    !
    !  Modified:
    !
    !    12 January 2008
    !
    !  Author:
    !
    !    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    KL Majumder, GP Bhattacharjee,
    !    Algorithm AS 63:
    !    The incomplete Beta Integral,
    !    Applied Statistics,
    !    Volume 22, Number 3, 1973, pages 409-411.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument, between 0 and 1.
    !
    !    Input, real ( kind = 8 ) P, Q, the parameters, which
    !    must be positive.
    !
    !    Input, real ( kind = 8 ) BETA, the logarithm of the complete
    !    beta function.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error flag.
    !    0, no error.
    !    nonzero, an error occurred.
    !
    !    Output, real ( kind = 8 ) BETAIN, the value of the incomplete
    !    Beta function ratio.
    !
      implicit none

      real ( kind = 8 ), parameter :: acu = 0.1D-14
      real ( kind = 8 ) ai
      real ( kind = 8 ) beta
      real ( kind = 8 ) betain
      real ( kind = 8 ) cx
      integer ( kind = 4 ) ifault
      logical indx
      integer ( kind = 4 ) ns
      real ( kind = 8 ) p
      real ( kind = 8 ) pp
      real ( kind = 8 ) psq
      real ( kind = 8 ) q
      real ( kind = 8 ) qq
      real ( kind = 8 ) rx
      real ( kind = 8 ) temp
      real ( kind = 8 ) term
      real ( kind = 8 ) x
      real ( kind = 8 ) xx

      betain = x
      ifault = 0
    !
    !  Check the input arguments.
    !
      if ( p <= 0.0D+00 .or. q <= 0.0D+00 ) then
        ifault = 1
        return
      end if

      if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
        ifault = 2
        return
      end if
    !
    !  Special cases.
    !
      if ( x == 0.0D+00 .or. x == 1.0D+00 ) then
        return
      end if
    !
    !  Change tail if necessary and determine S.
    !
      psq = p + q
      cx = 1.0D+00 - x

      if ( p < psq * x ) then
        xx = cx
        cx = x
        pp = q
        qq = p
        indx = .true.
      else
        xx = x
        pp = p
        qq = q
        indx = .false.
      end if

      term = 1.0D+00
      ai = 1.0D+00
      betain = 1.0D+00
      ns = int ( qq + cx * psq )
    !
    !  Use Soper's reduction formula.
    !
      rx = xx / cx
      temp = qq - ai
      if ( ns == 0 ) then
        rx = xx
      end if

      do

        term = term * temp * rx / ( pp + ai )
        betain = betain + term
        temp = abs ( term )

        if ( temp <= acu .and. temp <= acu * betain ) then

          betain = betain * exp ( pp * log ( xx ) &
            + ( qq - 1.0D+00 ) * log ( cx ) - beta ) / pp

          if ( indx ) then
            betain = 1.0D+00 - betain
          end if

          exit

        end if

        ai = ai + 1.0D+00
        ns = ns - 1

        if ( 0 <= ns ) then
          temp = qq - ai
          if ( ns == 0 ) then
            rx = xx
          end if
        else
          temp = psq
          psq = psq + 1.0D+00
        end if

      end do

      return
    end
    subroutine timestamp ( )

    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 May 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
      implicit none

      character ( len = 8 ) ampm
      integer ( kind = 4 ) d
      integer ( kind = 4 ) h
      integer ( kind = 4 ) m
      integer ( kind = 4 ) mm
      character ( len = 9 ), parameter, dimension(12) :: month = (/ &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' /)
      integer ( kind = 4 ) n
      integer ( kind = 4 ) s
      integer ( kind = 4 ) values(8)
      integer ( kind = 4 ) y

      call date_and_time ( values = values )

      y = values(1)
      m = values(2)
      d = values(3)
      h = values(5)
      n = values(6)
      s = values(7)
      mm = values(8)

      if ( h < 12 ) then
        ampm = 'AM'
      else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h < 12 ) then
          ampm = 'PM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
        d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

      return
    end

END MODULE MOD_IncBeta