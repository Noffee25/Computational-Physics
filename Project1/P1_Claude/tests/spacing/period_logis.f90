! =====================================================================
!  period_logis.f90 : eventual periodicity of the logistic map
!  x -> 4 x (1-x) iterated in IEEE double precision, exactly as in
!  rnglib (st = 4.0d0*st*(1.0d0-st)). Brent's cycle detection compares
!  the 64-bit patterns of the state.
!
!  Starts: the rnglib state for seed 777, followed by NSTART further
!  random starting points in (0,1) taken from MT19937 (53-bit).
!
!  The implementation is first checked against rnglib's outputs
!  (2/pi) asin(sqrt(x)) bit for bit for 10^6 steps.
!
!  usage: period_logis nstart mtseed outfile
!  output: start_index x0_bits x0 mu lambda min_state_on_cycle seconds
! =====================================================================
program period_logis
  use rnglib
  implicit none
  integer :: nstart, i, k
  integer(8) :: mtseed, lam, mu, t0, t1, rate
  real(8) :: x0, x, u, umod, cmin
  character(len=256) :: arg, outfile

  call get_command_argument(1, arg); read(arg, *) nstart
  call get_command_argument(2, arg); read(arg, *) mtseed
  call get_command_argument(3, outfile)

  ! ---- check against rnglib
  call rng_init(G_LOGIS, 777_8)
  x = modulo(dble(777_8) * 0.6180339887498949d0, 1.0d0)
  do k = 1, 1000000
     u = rng_next(G_LOGIS)
     x = 4.0d0 * x * (1.0d0 - x)
     umod = (2.0d0 / 3.141592653589793d0) * asin(sqrt(x))
     if (umod >= 1.0d0) umod = 1.0d0 - epsilon(1.0d0)
     if (transfer(umod, 1_8) /= transfer(u, 1_8)) then
        print *, 'mismatch at step', k, u, umod
        stop 1
     end if
  end do
  print '(a)', '# logistic map verified bitwise against rnglib for 1e6 steps'

  open(20, file=outfile, status='replace')
  write(20, '(a)') '# idx x0_bits x0 mu lambda cycle_min seconds'
  call rng_init(G_MT, mtseed)
  call system_clock(t0, rate)
  do i = 0, nstart
     if (i == 0) then
        x0 = modulo(dble(777_8) * 0.6180339887498949d0, 1.0d0)
     else
        x0 = 0.0d0
        do while (x0 <= 0.0d0)
           x0 = rng_next(G_MT)
        end do
     end if
     call system_clock(t0)
     call brent(x0, lam, mu, cmin)
     call system_clock(t1)
     write(20, '(i0,1x,i0,1x,es25.17,1x,i0,1x,i0,1x,es25.17,1x,f9.2)') i, transfer(x0, 1_8), x0, &
          mu, lam, cmin, dble(t1 - t0) / dble(rate)
     flush(20)
  end do
  close(20)

contains

  pure function f(x) result(y)
    real(8), intent(in) :: x
    real(8) :: y
    y = 4.0d0 * x * (1.0d0 - x)
  end function f

  subroutine brent(x0, lam, mu, cmin)
    real(8), intent(in) :: x0
    integer(8), intent(out) :: lam, mu
    real(8), intent(out) :: cmin
    integer(8) :: power, i
    real(8) :: tort, hare
    power = 1; lam = 1
    tort = x0; hare = f(x0)
    do while (transfer(tort, 1_8) /= transfer(hare, 1_8))
       if (power == lam) then
          tort = hare; power = 2 * power; lam = 0
       end if
       hare = f(hare); lam = lam + 1
    end do
    tort = x0; hare = x0
    do i = 1, lam
       hare = f(hare)
    end do
    mu = 0
    do while (transfer(tort, 1_8) /= transfer(hare, 1_8))
       tort = f(tort); hare = f(hare); mu = mu + 1
    end do
    ! smallest state on the cycle (canonical label of the cycle)
    cmin = tort
    hare = f(tort)
    do i = 1, lam - 1
       if (hare < cmin) cmin = hare
       hare = f(hare)
    end do
  end subroutine brent

end program period_logis
