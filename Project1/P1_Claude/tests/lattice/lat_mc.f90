! =====================================================================
! lat_mc.f90 : Monte Carlo integrals that the LCGs get exactly wrong.
!   I_s = integral over [0,1]^t of cos(2 pi (s.x - phi)) dx = 0  (s /= 0)
! For a dual-lattice vector s of a generator, every t-tuple satisfies
! s.x = phi (mod 1), so the MC estimate is exactly 1 with zero variance.
! Each (generator, integrand) pair uses N = 1e7 non-overlapping t-tuples,
! seed 777, and records the running mean and standard error at
! log-spaced checkpoints.
! =====================================================================
program lat_mc
  use rnglib
  implicit none
  integer, parameter :: NG = 7, NK = 8, NCP = 25
  integer(8), parameter :: NMAX = 10000000_8
  character(len=8) :: gname(NG) = [character(len=8) :: 'RANDU', 'RANDM', 'MINSTD', 'MINSTD2', &
       'LCG69069', 'MT19937', 'RANNUM']
  integer :: kt(NK), ks(8, NK), ig, ik, i, icp, unit, u2
  real(8) :: phi(NK), x(8), ph, f, mean, m2, delta, sx, sumc, sums
  integer(8) :: n, cp(NCP), st4, st5
  real(8), parameter :: TWOPI = 6.283185307179586d0

  kt = [3, 2, 3, 6, 2, 8, 8, 6]
  ks = 0
  ks(1:3, 1) = [9, -6, 1]
  ks(1:2, 2) = [6539, -1]
  ks(1:3, 3) = [25, 115, 678]
  ks(1:6, 4) = [3, 6, -8, -14, 9, -4]
  ks(1:2, 5) = [16807, -1]
  ks(1:8, 6) = [5, 4, 8, 1, 7, -2, 0, 1]
  ks(1:8, 7) = [3, -4, 4, 4, 2, -1, 2, 4]
  ks(1:6, 8) = [3, -11, 4, -8, 4, -4]
  phi = [0d0, 0d0, 0.5d0, 0.5d0, 0d0, 0d0, 0d0, 1917898033d0 / 4294967296d0]
  do icp = 1, NCP
     cp(icp) = nint(10.0d0 ** (1.0d0 + 0.25d0 * (icp - 1)), 8)
  end do

  open(newunit=unit, file='mc_running.txt', status='replace')
  write(unit, '(a)') '# gen killer N mean stderr   (integrand cos(2pi(s.x-phi)), exact integral 0)'
  open(newunit=u2, file='mc_final.txt', status='replace')
  write(u2, '(a)') '# gen killer t N mean_cos stderr |mean exp(2 pi i s.x)| z=mean/(1/sqrt(2N))'
  do ig = 1, NG
     do ik = 1, NK
        call init(ig)
        mean = 0; m2 = 0; icp = 1; sumc = 0; sums = 0
        do n = 1, NMAX
           do i = 1, kt(ik)
              x(i) = draw(ig)
           end do
           sx = 0
           do i = 1, kt(ik)
              sx = sx + dble(ks(i, ik)) * x(i)
           end do
           ph = TWOPI * (sx - phi(ik))
           f = cos(ph)
           sumc = sumc + cos(TWOPI * sx); sums = sums + sin(TWOPI * sx)
           delta = f - mean
           mean = mean + delta / dble(n)
           m2 = m2 + delta * (f - mean)
           if (icp <= NCP) then
              if (n == cp(icp)) then
                 write(unit, '(a9,i3,i10,es16.8,es14.6)') gname(ig), ik, n, mean, &
                      sqrt(m2 / dble(max(n - 1, 1_8)) / dble(n))
                 icp = icp + 1
              end if
           end if
        end do
        write(u2, '(a9,i3,i3,i10,es16.8,es14.6,f12.8,f14.2)') gname(ig), ik, kt(ik), NMAX, mean, &
             sqrt(m2 / dble(NMAX - 1) / dble(NMAX)), sqrt(sumc**2 + sums**2) / dble(NMAX), &
             mean / (1.0d0 / sqrt(2.0d0 * NMAX))
     end do
     print '(a,a)', 'done ', gname(ig)
  end do
  close(unit); close(u2)

contains

  subroutine init(ig)
    integer, intent(in) :: ig
    select case (ig)
    case (1); call rng_init(G_RANDU, 777_8)
    case (2); call rng_init(G_RANDM, 777_8)
    case (3); call rng_init(G_MINSTD, 777_8)
    case (4); st4 = 777_8
    case (5); st5 = 777_8
    case (6); call rng_init(G_MT, 777_8)
    case (7); call rng_init(G_RANNUM, 777_8)
    end select
  end subroutine init

  function draw(ig) result(y)
    integer, intent(in) :: ig
    real(8) :: y
    select case (ig)
    case (1); y = rng_next(G_RANDU)
    case (2); y = rng_next(G_RANDM)
    case (3); y = rng_next(G_MINSTD)
    case (4)
       st4 = modulo(48271_8 * st4, 2147483647_8); y = dble(st4) / 2147483647d0
    case (5)
       st5 = iand(69069_8 * st5 + 1_8, 4294967295_8); y = dble(st5) / 4294967296d0
    case (6); y = rng_next(G_MT)
    case (7); y = rng_next(G_RANNUM)
    case default; y = 0
    end select
  end function draw

end program lat_mc
