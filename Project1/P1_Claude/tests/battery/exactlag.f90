! =====================================================================
!  exactlag.f90 : exact full-period serial correlations of RANDU/RANDM
!  using 128-bit integer sums over one complete period (2^29 states).
!  Because the sequence is periodic, the cyclic sum over one period is
!  the exact population value of the lag-k correlation:
!      r_k = (P sum s_n s_{n+k} - (sum s)^2) / (P sum s^2 - (sum s)^2)
!  with s_{n+k} = a^k s_n mod 2^31.  Also verifies the period and the
!  half-period identity s_{n+2^28} = s_n + 2^30 (mod 2^31).
!  usage: exactlag a seed
! =====================================================================
program exactlag
  implicit none
  integer, parameter :: NL = 6
  integer(8), parameter :: M31 = 2147483648_8, MSK = 2147483647_8
  integer :: lags(NL) = [1, 2, 3, 4, 65, 256]
  integer(8) :: a, seed, s, s0, b(NL), y, P, n, half
  integer(16) :: sx, sxx, sxy(NL), PP
  real(16) :: r, num, den
  integer :: i, k, nbad
  character(len=64) :: arg
  call get_command_argument(1, arg); read(arg,*) a
  call get_command_argument(2, arg); read(arg,*) seed
  do i = 1, NL
     b(i) = 1
     do k = 1, lags(i); b(i) = iand(b(i) * a, MSK); end do
  end do
  ! half-period multiplier a^(2^28)
  half = a
  do k = 1, 28; half = iand(half*half, MSK); end do
  print '(a,i12,a,i12)', 'a^(2^28) mod 2^31 = ', half, '   2^30+1 = ', 1073741825_8
  s0 = modulo(seed, M31); s = s0
  sx = 0; sxx = 0; sxy = 0
  P = 0; nbad = 0
  do
     sx = sx + s
     sxx = sxx + int(s,16)*s
     do i = 1, NL
        y = iand(b(i) * s, MSK)
        sxy(i) = sxy(i) + int(s,16)*y
     end do
     y = iand(half * s, MSK)
     if (y /= iand(s + 1073741824_8, MSK)) nbad = nbad + 1
     s = iand(a * s, MSK)
     P = P + 1
     if (s == s0) exit
  end do
  print '(a,i12,a,i12)', 'period = ', P, '   2^29 = ', 536870912_8
  print '(a,i12)', 'violations of s_{n+P/2} = s_n + 2^30 mod 2^31 : ', nbad
  PP = P
  print '(a,es40.30)', 'mean x       = ', real(sx,16)/real(PP,16)/real(M31,16)
  den = real(PP*sxx - sx*sx, 16)
  do i = 1, NL
     num = real(PP*sxy(i) - sx*sx, 16)
     r = num/den
     print '(a,i4,a,i12,a,es30.20,a,f9.4)', 'lag ', lags(i), '  a^k mod m = ', b(i), '  r_k = ', r, &
          '  sqrt(P)*r_k = ', real(r*sqrt(real(PP,16)), 8)
  end do
end program
