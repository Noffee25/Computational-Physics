! =====================================================================
! hist_period.f90 : "population" histogram of one full period
!
! For a generator that visits each of its reachable states exactly once
! per period P, the histogram of ONE FULL PERIOD gives c_b, the number
! of reachable values that land in bin b.  Every finite-sample histogram
! of that generator is then a (nearly) deterministic function of c_b:
! N = q P + r draws give n_b = q c_b + h_b with h_b ~ Hypergeometric(P,
! c_b, r).  This program bins NPER draws into B = 10^1 ... 10^8 bins
! simultaneously, using exactly the arithmetic of rnd_test_claude.f, and
! writes the count-of-counts table of c_b for each depth.
!
! It also verifies the period: for generators with an integer state the
! state must return to its initial value after exactly NPER steps and
! never before; for RAND the 16 outputs after NPER draws must repeat the
! first 16.
!
! usage: hist_period IGEN IPREC SEED NPER TAG OUTPREFIX
! =====================================================================
program hist_period
  use rnglib
  implicit none
  type cnt_t
     integer, allocatable :: c(:)
  end type cnt_t
  type(cnt_t) :: h(8)
  integer :: igen, iprec, d, nb(8), idx, k, mn, mx, nret
  integer(8) :: seed, nper, i, s0, tot
  integer(8), allocatable :: tab(:)
  real(8) :: x, db(8), first(16), y
  character(len=64) :: arg, tag, pre
  logical :: okrep
  integer :: t0, t1, rate

  call get_command_argument(1, arg); read(arg, *) igen
  call get_command_argument(2, arg); read(arg, *) iprec
  call get_command_argument(3, arg); read(arg, *) seed
  call get_command_argument(4, arg); read(arg, *) nper
  call get_command_argument(5, tag)
  call get_command_argument(6, pre)
  call system_clock(t0, rate)

  do d = 1, 8
     nb(d) = 10**d
     db(d) = dble(nb(d))
     allocate(h(d)%c(nb(d)))
     h(d)%c = 0
  end do

  rng_single = (iprec == 4)
  call rng_init(igen, seed)
  s0 = rng_state_int(igen)
  nret = 0
  do i = 1, nper
     x = rng_next(igen)
     if (i <= 16) first(i) = x
     do d = 1, 8
        idx = int(x * db(d)) + 1
        if (idx > nb(d)) idx = nb(d)
        if (idx < 1) idx = 1
        h(d)%c(idx) = h(d)%c(idx) + 1
     end do
     if (s0 >= 0) then
        if (rng_state_int(igen) == s0) nret = nret + 1
     end if
  end do
  ! period check
  okrep = .true.
  do i = 1, 16
     y = rng_next(igen)
     if (y /= first(i)) okrep = .false.
  end do
  rng_single = .false.
  call system_clock(t1)

  open(10, file=trim(pre)//'_coc.out', status='unknown', position='append')
  open(11, file=trim(pre)//'_sum.out', status='unknown', position='append')
  do d = 1, 8
     mn = minval(h(d)%c); mx = maxval(h(d)%c)
     tot = sum(int(h(d)%c, 8))
     if (mx - mn > 50000000) stop 'c range too wide'
     allocate(tab(mn:mx)); tab = 0
     do k = 1, nb(d)
        tab(h(d)%c(k)) = tab(h(d)%c(k)) + 1
     end do
     do k = mn, mx
        if (tab(k) > 0) write(10, '(A,1X,I1,1X,I11,1X,I11)') trim(tag), d, k, tab(k)
     end do
     write(11, '(A,1X,I2,1X,I1,1X,I12,1X,I10,1X,I12,1X,I12,2(1X,I11),1X,I6,1X,L1,1X,F8.1)') &
          trim(tag), igen, iprec, seed, nb(d), nper, tot, mn, mx, nret, okrep, &
          dble(t1 - t0) / dble(rate)
     deallocate(tab)
  end do
  close(10); close(11)
end program hist_period
