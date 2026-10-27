! (1) Checks, draw by draw, that gfortran RAND() after srand(seed) equals
!     floor((s-1)/2^9)/2^22 where s is the Park-Miller (MINSTD) state from
!     the same seed, and compares with the alternative floor(2^22 s/M)/2^22.
! (2) Builds the exact one-period population table c_b of RAND for
!     B = 10^1..10^8: each grid value k/2^22 carries 512 MINSTD states
!     (the top one 510), binned with the arithmetic of rnd_test_claude.f.
program rand_check
  use rnglib
  implicit none
  integer(8), parameter :: NCHK = 200000000_8, M = 2147483647_8
  integer(8) :: i, s, nbad1, nbad2, tot, k
  integer, allocatable :: c(:)
  integer(8), allocatable :: tab(:)
  real(8) :: x, y1, y2, db
  integer :: d, nb, idx, mn, mx, w, j
  call rng_init(G_RAND, 777_8); call rng_init(G_MINSTD, 777_8)
  nbad1 = 0; nbad2 = 0
  do i = 1, NCHK
     x = rng_next(G_RAND)
     y1 = rng_next(G_MINSTD)
     s = rng_state_int(G_MINSTD)
     y1 = dble((s - 1) / 512_8) / 4194304.0d0
     y2 = dble((s * 4194304_8) / M) / 4194304.0d0
     if (x /= y1) nbad1 = nbad1 + 1
     if (x /= y2) nbad2 = nbad2 + 1
  end do
  write(*,'(A,I12,A,I10,A,I10)') 'draws checked', NCHK, '  mismatches floor((s-1)/512):', nbad1, &
       '  mismatches floor(2^22 s/M):', nbad2
  open(10, file='per_coc.out', status='old', position='append')
  open(11, file='per_sum.out', status='old', position='append')
  do d = 1, 8
     nb = 10**d; db = dble(nb)
     allocate(c(nb)); c = 0
     do k = 0, 4194303_8
        w = 512; if (k == 4194303_8) w = 510
        x = dble(k) / 4194304.0d0
        idx = int(x * db) + 1
        if (idx > nb) idx = nb
        c(idx) = c(idx) + w
     end do
     mn = minval(c); mx = maxval(c); tot = sum(int(c, 8))
     allocate(tab(mn:mx)); tab = 0
     do j = 1, nb
        tab(c(j)) = tab(c(j)) + 1
     end do
     do j = mn, mx
        if (tab(j) > 0) write(10, '(A,1X,I1,1X,I11,1X,I11)') 'RAND_d', d, j, tab(j)
     end do
     write(11, '(A,1X,I2,1X,I1,1X,I12,1X,I10,1X,I12,1X,I12,2(1X,I11),1X,I6,1X,L1,1X,F8.1)') &
          'RAND_d', 3, 8, 777_8, nb, 2147483646_8, tot, mn, mx, -1, .true., 0.0
     deallocate(c, tab)
  end do
end program rand_check
