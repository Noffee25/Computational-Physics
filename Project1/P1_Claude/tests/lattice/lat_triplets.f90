! =====================================================================
! lat_triplets.f90 : RANDU triplets and the 15 planes 9x - 6y + z = j.
!  * 1e7 overlapping triplets (x_k,x_{k+1},x_{k+2}) from seed 777:
!    exact integer plane index j = (9X_k - 6X_{k+1} + X_{k+2}) / 2^31,
!    plane populations, check that the division is always exact
!  * 1e7 non-overlapping triplets, same statistics
!  * floating point view: S = 9x - 6y + z computed from the REAL(8) outputs,
!    and with rng_single (outputs rounded to REAL*4)
!  * the same statistic for MT19937 and RANNUM (histogram on [-6,10])
!  * small samples of triplets for scatter plots
! =====================================================================
program lat_triplets
  use rnglib
  implicit none
  integer(8), parameter :: N = 10000000_8, M = 2147483648_8
  integer, parameter :: NB = 1600            ! bins of width 0.01 on [-6,10)
  integer(8) :: X(3), r, j, cnt_ov(-6:10), cnt_no(-6:10), nbad, k
  integer(8) :: hist(NB, 4), nnear(4)
  real(8) :: x3(3), sv, dev, devmax(2), devhist(0:40)
  integer :: ig, ids(4), ib, i, unit
  character(len=8) :: nm(4)

  ! ---------------- exact integer planes, overlapping triplets
  call rng_init(G_RANDU, 777_8)
  sv = rng_next(G_RANDU); X(1) = rng_state_int(G_RANDU)
  sv = rng_next(G_RANDU); X(2) = rng_state_int(G_RANDU)
  cnt_ov = 0; nbad = 0
  do k = 1, N
     sv = rng_next(G_RANDU); X(3) = rng_state_int(G_RANDU)
     r = 9 * X(1) - 6 * X(2) + X(3)
     if (modulo(r, M) /= 0) nbad = nbad + 1
     j = r / M
     cnt_ov(j) = cnt_ov(j) + 1
     X(1) = X(2); X(2) = X(3)
  end do
  print '(a,i0)', 'overlapping triplets with 9X-6Y+Z not divisible by 2^31: ', nbad
  ! ---------------- non-overlapping triplets
  call rng_init(G_RANDU, 777_8)
  cnt_no = 0; nbad = 0
  do k = 1, N
     do i = 1, 3
        sv = rng_next(G_RANDU); X(i) = rng_state_int(G_RANDU)
     end do
     r = 9 * X(1) - 6 * X(2) + X(3)
     if (modulo(r, M) /= 0) nbad = nbad + 1
     cnt_no(r / M) = cnt_no(r / M) + 1
  end do
  print '(a,i0)', 'non-overlapping triplets not divisible: ', nbad
  open(newunit=unit, file='triplet_planes.txt', status='replace')
  write(unit, '(a)') '# j  count_overlapping  count_nonoverlapping   (RANDU seed 777, 1e7 triplets each)'
  do j = -6, 10
     write(unit, '(i4,2i12)') j, cnt_ov(j), cnt_no(j)
  end do
  close(unit)

  ! ---------------- floating-point view and the other generators
  ids = [G_RANDU, G_MT, G_RANNUM, G_RANDU]
  nm  = [character(len=8) :: 'RANDU', 'MT19937', 'RANNUM', 'RANDU_SP']
  hist = 0; nnear = 0; devmax = 0; devhist = 0
  do ig = 1, 4
     rng_single = (ig == 4)
     call rng_init(ids(ig), 777_8)
     do k = 1, N
        x3(1) = rng_next(ids(ig)); x3(2) = rng_next(ids(ig)); x3(3) = rng_next(ids(ig))
        sv = 9.0d0 * x3(1) - 6.0d0 * x3(2) + x3(3)
        dev = abs(sv - anint(sv))
        if (dev < 0.005d0) nnear(ig) = nnear(ig) + 1
        if (ig == 1) devmax(1) = max(devmax(1), dev)
        if (ig == 4) then
           devmax(2) = max(devmax(2), dev)
           ! log10 histogram of the rounding residual
           if (dev > 0) then
              ib = int(-log10(dev) * 2.0d0)
              ib = max(0, min(40, ib))
           else
              ib = 40
           end if
           devhist(ib) = devhist(ib) + 1
        end if
        ib = int((sv + 6.0d0 + 0.005d0) / 0.01d0) + 1   ! bins centred on multiples of 0.01
        if (ib >= 1 .and. ib <= NB) hist(ib, ig) = hist(ib, ig) + 1
     end do
  end do
  rng_single = .false.
  print '(a,es12.4)', 'RANDU double: max |S - nint(S)| = ', devmax(1)
  print '(a,es12.4)', 'RANDU single: max |S - nint(S)| = ', devmax(2)
  do ig = 1, 4
     print '(a8,a,i10,a,f9.6)', nm(ig), ' triplets with |S-nint(S)|<0.005: ', nnear(ig), '  fraction ', dble(nnear(ig)) / N
  end do
  open(newunit=unit, file='triplet_hist.txt', status='replace')
  write(unit, '(a)') '# bin_center  RANDU  MT19937  RANNUM  RANDU_SP  (1e7 non-overlapping triplets, S = 9x-6y+z)'
  do ib = 1, NB
     write(unit, '(f9.3,4i10)') -6.0d0 + 0.01d0 * (ib - 1), hist(ib, :)
  end do
  close(unit)
  open(newunit=unit, file='triplet_sp_resid.txt', status='replace')
  write(unit, '(a)') '# k  count of single-precision residual |S-nint(S)| in [10^-(k+1)/2, 10^-k/2)'
  do ib = 0, 40
     write(unit, '(i4,f14.1)') ib, devhist(ib)
  end do
  close(unit)

  ! ---------------- scatter samples
  open(newunit=unit, file='triplet_sample.txt', status='replace')
  write(unit, '(a)') '# RANDU x y z   MT19937 x y z   (30000 non-overlapping triplets, seed 777)'
  call rng_init(G_RANDU, 777_8); call rng_init(G_MT, 777_8)
  do k = 1, 30000
     x3(1) = rng_next(G_RANDU); x3(2) = rng_next(G_RANDU); x3(3) = rng_next(G_RANDU)
     write(unit, '(3f11.7)', advance='no') x3
     x3(1) = rng_next(G_MT); x3(2) = rng_next(G_MT); x3(3) = rng_next(G_MT)
     write(unit, '(3f11.7)') x3
  end do
  close(unit)
end program lat_triplets
