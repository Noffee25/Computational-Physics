! =====================================================================
!  knuth.f90 : short battery of classical empirical tests (Knuth 3.3.2)
!  Each test consumes its own consecutive segment of the stream and
!  writes a table of observed counts; probabilities and chi-square
!  p-values are computed afterwards (battery_analysis.py).
!    RUNS : runs up, element after each run discarded (independent runs)
!           lengths 1..5, >=6                      (Knuth 3.3.2 G, ex. 14)
!    GAP  : gap test, [alpha,beta) = [0,0.1), gaps 0..39, >=40 (3.3.2 B)
!    PERM : permutation test, non-overlapping 4-tuples, 24 orders (3.3.2 F)
!    COUP : coupon collector, d = 8, lengths 8..39, >=40   (3.3.2 D)
!    MAXT : maximum of t = 8, V = max^8 in 100 equal bins  (3.3.2 H)
!    BDAY : birthday spacings, m = 512 birthdays in n = 2^24 days,
!           R = 16384 samples, J = number of repeated spacings (3.3.2 J)
!    LAG1 : sum of P1(x_n)P1(x_{n+1}) over the segment (serial corr.)
!  variant 0: u = x ; variant 1 (low bits): u = frac(2^8 x), i.e. the
!  leading 8 bits are discarded
!  usage: knuth gen seed variant log2n outfile
! =====================================================================
program knuth
  use rnglib
  implicit none
  integer :: gen, variant, lg, i, r, len, j, q, ntie
  integer(8) :: seed, n, cnt, used
  integer(8) :: runs(6), gaps(0:40), perm(0:23), coup(8:40), maxt(100), bday(0:7)
  integer :: days(512), sp(512), nrep, ncoup_fail
  real(8) :: u, prev, v(8), mx, s, plag, pprev
  logical :: seen(0:7)
  integer :: nseen
  character(len=256) :: arg, outfile
  call get_command_argument(1, arg); read(arg,*) gen
  call get_command_argument(2, arg); read(arg,*) seed
  call get_command_argument(3, arg); read(arg,*) variant
  call get_command_argument(4, arg); read(arg,*) lg
  call get_command_argument(5, outfile)
  n = 2_8**lg
  call rng_init(gen, seed)
  open(11, file=trim(outfile), status='replace')
  write(11,'(a,a,a,i12,a,i2,a,i14)') '# ', gen_name(gen), ' seed ', seed, ' variant ', variant, ' n ', n

  ! ---------------- RUNS UP (independent) ----------------
  runs = 0; used = 0; ntie = 0
  do while (used < n)
     prev = nxt(); len = 1
     do
        u = nxt()
        if (u > prev) then
           len = len + 1; prev = u
        else
           if (u == prev) ntie = ntie + 1
           exit                         ! u (the run breaker) is discarded
        end if
     end do
     runs(min(len,6)) = runs(min(len,6)) + 1
  end do
  write(11,'(a,i10)') 'RUNS ties ', ntie
  write(11,'(6i12)') runs

  ! ---------------- GAP [0,0.1) ----------------
  gaps = 0; used = 0; r = 0
  do while (used < n)
     u = nxt()
     if (u < 0.1d0) then
        gaps(min(r,40)) = gaps(min(r,40)) + 1; r = 0
     else
        r = r + 1
     end if
  end do
  write(11,'(a)') 'GAP'
  write(11,'(41i12)') gaps

  ! ---------------- PERMUTATION t=4 ----------------
  perm = 0; used = 0
  do while (used + 4 <= n)
     do i = 1, 4; v(i) = nxt(); end do
     q = lehmer4(v)
     perm(q) = perm(q) + 1
  end do
  write(11,'(a)') 'PERM'
  write(11,'(24i12)') perm

  ! ---------------- COUPON COLLECTOR d=8 ----------------
  coup = 0; used = 0; ncoup_fail = 0
  outer: do while (used < n)
     seen = .false.; nseen = 0; len = 0
     do while (nseen < 8)
        u = nxt(); len = len + 1
        q = min(int(8.0d0*u), 7)
        if (.not. seen(q)) then
           seen(q) = .true.; nseen = nseen + 1
        end if
        if (len > 100000) then          ! never completes (degenerate generator)
           ncoup_fail = ncoup_fail + 1
           if (used >= n) exit outer
           len = 0; seen = .false.; nseen = 0
        end if
        if (used >= n .and. nseen < 8) exit outer
     end do
     coup(min(len,40)) = coup(min(len,40)) + 1
  end do outer
  write(11,'(a,i10)') 'COUP fails ', ncoup_fail
  write(11,'(33i12)') coup

  ! ---------------- MAXIMUM OF 8 ----------------
  maxt = 0; used = 0
  do while (used + 8 <= n)
     mx = 0
     do i = 1, 8; u = nxt(); mx = max(mx, u); end do
     q = min(int(100.0d0 * mx**8), 99) + 1
     maxt(q) = maxt(q) + 1
  end do
  write(11,'(a)') 'MAXT'
  write(11,'(100i10)') maxt

  ! ---------------- BIRTHDAY SPACINGS ----------------
  bday = 0
  do j = 1, 16384
     do i = 1, 512
        u = nxt()
        days(i) = min(int(u * 16777216.0d0), 16777215)
     end do
     call isort(days, 512)
     do i = 2, 512; sp(i) = days(i) - days(i-1); end do
     sp(1) = days(1) + 16777216 - days(512)
     call isort(sp, 512)
     nrep = 0
     do i = 2, 512
        if (sp(i) == sp(i-1)) nrep = nrep + 1
     end do
     bday(min(nrep,7)) = bday(min(nrep,7)) + 1
  end do
  write(11,'(a)') 'BDAY'
  write(11,'(8i12)') bday

  ! ---------------- LAG-1 SERIAL ----------------
  used = 0; s = 0
  pprev = 1.7320508075688772d0 * (2.0d0*nxt() - 1.0d0)
  cnt = 0
  do while (used < n)
     plag = 1.7320508075688772d0 * (2.0d0*nxt() - 1.0d0)
     s = s + plag * pprev; pprev = plag; cnt = cnt + 1
  end do
  write(11,'(a)') 'LAG1'
  write(11,'(es24.15,i14)') s, cnt
  close(11)

contains
  real(8) function nxt()
    real(8) :: xx
    xx = rng_next(gen)
    if (variant == 1) then
       xx = xx * 256.0d0
       xx = xx - aint(xx)
    end if
    used = used + 1
    nxt = xx
  end function

  integer function lehmer4(w)
    ! index 0..23 of the relative order of w(1..4); ties broken by position
    real(8), intent(in) :: w(8)
    integer :: a, b, c
    lehmer4 = 0
    do a = 1, 4
       c = 0
       do b = a+1, 4
          if (w(b) < w(a)) c = c + 1
       end do
       lehmer4 = lehmer4 * (5 - a) + c
    end do
  end function

  subroutine isort(a, m)
    ! heapsort of integer array
    integer, intent(inout) :: a(:)
    integer, intent(in) :: m
    integer :: ii, jj, ll, ir, ra
    if (m < 2) return
    ll = m/2 + 1; ir = m
    do
       if (ll > 1) then
          ll = ll - 1; ra = a(ll)
       else
          ra = a(ir); a(ir) = a(1); ir = ir - 1
          if (ir == 1) then
             a(1) = ra; return
          end if
       end if
       ii = ll; jj = ll + ll
       do while (jj <= ir)
          if (jj < ir) then
             if (a(jj) < a(jj+1)) jj = jj + 1
          end if
          if (ra < a(jj)) then
             a(ii) = a(jj); ii = jj; jj = jj + jj
          else
             jj = ir + 1
          end if
       end do
       a(ii) = ra
    end do
  end subroutine
end program
