program logcycle
  ! Brent cycle detection on the double-precision logistic map state x -> 4x(1-x)
  implicit none
  real(8) :: tort, hare
  integer(8) :: power, lam, mu, i, maxit
  real(8) :: x0
  x0 = modulo(777d0 * 0.6180339887498949d0, 1.0d0)
  maxit = 3000000000_8
  power = 1; lam = 1
  tort = x0; hare = f(x0)
  do while (tort /= hare)
     if (power == lam) then
        tort = hare; power = power * 2; lam = 0
     end if
     hare = f(hare); lam = lam + 1
     if (power > maxit) then
        print *, 'no cycle found up to', maxit; stop
     end if
  end do
  print *, 'cycle length', lam
  ! find mu
  tort = x0; hare = x0
  do i = 1, lam; hare = f(hare); end do
  mu = 0
  do while (tort /= hare)
     tort = f(tort); hare = f(hare); mu = mu + 1
  end do
  print *, 'transient length', mu, ' cycle value', tort
contains
  real(8) function f(x)
    real(8), intent(in) :: x
    f = 4.0d0 * x * (1.0d0 - x)
  end function
end program
