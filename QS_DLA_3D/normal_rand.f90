module normal_rand
implicit none

contains


function r8_uniform_01 ( seed )

  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function r8_uniform_01

  function r8_normal_ab ( a, b, seed )

!*****************************************************************************80
!
!! R8_NORMAL_AB returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_ab
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_ab = a + b * x

  return
  end function r8_normal_ab



end program