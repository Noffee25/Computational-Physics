! lag-1 correlation of bit 49 of the MINSTD double, of bit 18, and of the
! rounding-carry flip c = b49 XOR b18 (copy violation), N draws
program minstd_carry
  use rnglib
  implicit none
  integer(8) :: n, nd, m, seed
  integer :: b49, b18, c, p49, p18, pc
  real(8) :: x, s(3), s11(3), x1, x2
  character(len=64) :: arg
  call get_command_argument(1, arg); read(arg, *) nd
  call get_command_argument(2, arg); read(arg, *) seed
  call rng_init(G_MINSTD, seed)
  s = 0; s11 = 0; p49 = 0; p18 = 0; pc = 0; x1 = 0; x2 = 0
  do n = 1, nd
     x = rng_next(G_MINSTD)
     m = int(x * 9007199254740992.0d0, 8)
     b49 = int(ibits(m, 53 - 49, 1)); b18 = int(ibits(m, 53 - 18, 1)); c = ieor(b49, b18)
     s(1) = s(1) + b49; s(2) = s(2) + b18; s(3) = s(3) + c
     if (n > 1) then
        s11(1) = s11(1) + b49 * p49; s11(2) = s11(2) + b18 * p18; s11(3) = s11(3) + c * pc
        x1 = x1 + pc * b18; x2 = x2 + p18 * c
     end if
     p49 = b49; p18 = b18; pc = c
  end do
  s = s / nd; s11 = s11 / (nd - 1)
  print '(a,3f12.6)', 'mean b49 b18 c   ', s
  print '(a,3f12.6)', 'lag1 corr        ', (s11 - s * s) / (s * (1 - s))
  print '(a,3f12.2)', 'z = r sqrt(N)    ', (s11 - s * s) / (s * (1 - s)) * sqrt(dble(nd))
  x1 = x1 / (nd - 1); x2 = x2 / (nd - 1)
  print '(a,2f12.2)', 'z cross c(n)b18(n+1), b18(n)c(n+1)', &
       (x1 - s(3) * s(2)) / sqrt(s(3) * (1 - s(3)) * s(2) * (1 - s(2))) * sqrt(dble(nd)), &
       (x2 - s(3) * s(2)) / sqrt(s(3) * (1 - s(3)) * s(2) * (1 - s(2))) * sqrt(dble(nd))
end program minstd_carry
