! Brent cycle detection for the floating point logistic map used by rnglib
program logcycle
  use rnglib
  implicit none
  real(8) :: s0, tort, hare, s, x, y
  integer(8) :: power, lam, mu, n, nmax, i
  call rng_init(G_LOGIS, 777_8)
  s0 = modulo(dble(777_8) * 0.6180339887498949d0, 1.0d0)
  ! check the replica against rnglib
  s = s0
  do i = 1, 1000000
     s = 4.0d0 * s * (1.0d0 - s)
     y = (2.0d0 / 3.141592653589793d0) * asin(sqrt(s))
     x = rng_next(G_LOGIS)
     if (x /= y) then
        print *, 'replica mismatch at', i; stop
     end if
  end do
  print *, 'replica matches rnglib for 1e6 steps'
  nmax = 4000000000_8
  power = 1; lam = 1
  tort = s0; hare = f(s0); n = 1
  do while (tort /= hare)
     if (power == lam) then
        tort = hare; power = power * 2; lam = 0
     end if
     hare = f(hare); lam = lam + 1; n = n + 1
     if (n > nmax) then
        print *, 'no cycle within', nmax, ' steps'; stop
     end if
  end do
  tort = s0; hare = s0
  do i = 1, lam
     hare = f(hare)
  end do
  mu = 0
  do while (tort /= hare)
     tort = f(tort); hare = f(hare); mu = mu + 1
  end do
  print *, 'cycle length', lam, ' tail', mu, ' value', tort
contains
  real(8) function f(z)
    real(8), intent(in) :: z
    f = 4.0d0 * z * (1.0d0 - z)
  end function f
end program logcycle
