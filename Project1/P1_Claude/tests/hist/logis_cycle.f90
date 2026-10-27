! Brent cycle detection for the double-precision logistic map used by
! rnglib's G_LOGIS (state s -> 4 s (1-s)).  First verifies that a local
! replica reproduces rng_next(G_LOGIS) bit for bit, then finds the tail
! length mu and the cycle length lam of the state orbit from seed 777
! (and a few other seeds given on the command line).
program logis_cycle
  use rnglib
  implicit none
  integer(8) :: seed, i, lam, mu, power, nmax
  real(8) :: s, t, h, x, y
  integer :: iarg, nbad
  character(len=32) :: arg
  nmax = 4000000000_8
  do iarg = 1, max(1, command_argument_count())
     seed = 777_8
     if (command_argument_count() >= 1) then
        call get_command_argument(iarg, arg); read(arg, *) seed
     end if
     s = st0(seed)
     ! --- verify replica against the library for 10^7 draws
     call rng_init(G_LOGIS, seed)
     nbad = 0; t = s
     do i = 1, 10000000_8
        t = f(t)
        x = (2.0d0 / 3.141592653589793d0) * asin(sqrt(t))
        if (x >= 1.0d0) x = 1.0d0 - epsilon(1.0d0)
        y = rng_next(G_LOGIS)
        if (x /= y) nbad = nbad + 1
     end do
     ! --- Brent: find lam
     power = 1; lam = 1
     t = s; h = f(s)
     do while (t /= h)
        if (power == lam) then
           t = h; power = power * 2; lam = 0
        end if
        h = f(h); lam = lam + 1
        if (lam > nmax) exit
     end do
     ! --- find mu
     t = s; h = s
     do i = 1, lam
        h = f(h)
     end do
     mu = 0
     do while (t /= h)
        t = f(t); h = f(h); mu = mu + 1
     end do
     write(*,'(A,I12,A,I8,A,I14,A,I14,A,ES24.16)') 'seed', seed, ' mismatches', nbad, &
          '  tail mu', mu, '  cycle lam', lam, '  cycle entry state', t
  end do
contains
  real(8) function f(z)
    real(8), intent(in) :: z
    f = 4.0d0 * z * (1.0d0 - z)
  end function f
  real(8) function st0(sd)
    integer(8), intent(in) :: sd
    st0 = modulo(dble(sd) * 0.6180339887498949d0, 1.0d0)
    if (st0 <= 0.0d0 .or. st0 >= 1.0d0) st0 = 0.3141592653589793d0
  end function st0
end program logis_cycle
