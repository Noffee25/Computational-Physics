! Hit-or-miss estimate of pi from the volume of the unit ball:
! v = 2u - 1 from non-overlapping triplets, pi_hat = 6 * (fraction with |v| <= 1).
program lat_pi
  use rnglib
  implicit none
  integer, parameter :: NG = 5, NCP = 33
  integer :: ids(NG), ig, icp, unit
  integer(8) :: n, nin, cp(NCP), nmax
  real(8) :: v(3), p
  character(len=8) :: nm(NG)
  ids = [G_RANDU, G_RANDM, G_MINSTD, G_MT, G_RANNUM]
  nm = [character(len=8) :: 'RANDU', 'RANDM', 'MINSTD', 'MT19937', 'RANNUM']
  do icp = 1, NCP
     cp(icp) = nint(10d0 ** (1d0 + 0.25d0 * (icp - 1)), 8)
  end do
  nmax = cp(NCP)   ! 1e9 is too much: stop at 2e8
  nmax = 200000000_8
  open(newunit=unit, file='pi_running.txt', status='replace')
  write(unit, '(a)') '# gen N pi_hat stderr   (pi_hat = 6 * fraction of triplets inside the unit ball)'
  do ig = 1, NG
     call rng_init(ids(ig), 777_8)
     nin = 0; icp = 1
     do n = 1, nmax
        v(1) = 2d0 * rng_next(ids(ig)) - 1d0
        v(2) = 2d0 * rng_next(ids(ig)) - 1d0
        v(3) = 2d0 * rng_next(ids(ig)) - 1d0
        if (v(1)**2 + v(2)**2 + v(3)**2 <= 1d0) nin = nin + 1
        if (icp <= NCP) then
           if (n == cp(icp) .or. n == nmax) then
              p = dble(nin) / dble(n)
              write(unit, '(a9,i11,f14.8,es13.5)') nm(ig), n, 6d0 * p, 6d0 * sqrt(p * (1 - p) / dble(n))
              if (n == cp(icp)) icp = icp + 1
           end if
        end if
     end do
     print '(a8,a,f12.7,a,es10.3)', nm(ig), ' pi_hat(2e8) = ', 6d0 * dble(nin) / dble(nmax), &
          '  stderr ', 6d0 * sqrt(p * (1 - p) / dble(nmax))
  end do
  close(unit)
end program
