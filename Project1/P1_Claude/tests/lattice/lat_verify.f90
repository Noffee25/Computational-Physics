! =====================================================================
! lat_verify.f90 : exhaustive, exact (integer) check of the spectral-test
! predictions.  For one generator it runs over the FULL period, forms every
! overlapping t-tuple (X_n,...,X_{n+t-1}) of integer states, and for each
! predicted normal vector s (read from vectors.txt) computes r = s.X exactly.
! It records which residues r mod M occur, how many distinct hyperplanes
! (distinct values of r) are occupied, the population of each plane and the
! largest gap between adjacent occupied planes, in unit-cube distance
! gap / (M |s|).
!
! usage: lat_verify <gen> [nsteps]   gen: 1 RANDU 2 RANDM 3 MINSTD
!                                         4 MINSTD2(48271) 5 LCG69069(c=1)
! =====================================================================
program lat_verify
  use rnglib
  implicit none
  integer, parameter :: MAXV = 32, MAXRES = 8
  integer :: gen, nv, iv, t, i, ios, ng, mode(MAXV), tv(MAXV)
  integer(8) :: s(8, MAXV), hnum(MAXV), hden(MAXV), gran(MAXV)
  integer(8) :: M, a, c, x, x0, period, nsteps, n, r, j, rr, win(8)
  integer(8) :: jlo(MAXV), jhi(MAXV), nres(MAXV), res(MAXRES, MAXV), l1
  integer(8), allocatable :: cnt(:,:,:)
  logical :: overflow(MAXV), back
  character(len=200) :: line, arg
  character(len=12) :: gname
  real(8) :: snorm, gapmax, gapu, t0, t1
  integer(8), allocatable :: pos(:)
  integer :: np, k, kk
  integer(8) :: lo2, ir, x7
  logical :: halfret

  call get_command_argument(1, arg); read(arg, *) gen
  nsteps = -1
  if (command_argument_count() >= 2) then
     call get_command_argument(2, arg); read(arg, *) nsteps
  end if
  x0 = 777_8
  select case (gen)
  case (1); M = 2147483648_8; a = 65539_8; c = 0; period = 536870912_8; gname = 'RANDU'
  case (2); M = 2147483648_8; a = 6539_8;  c = 0; period = 536870912_8; gname = 'RANDM'
  case (3); M = 2147483647_8; a = 16807_8; c = 0; period = 2147483646_8; gname = 'MINSTD'
  case (4); M = 2147483647_8; a = 48271_8; c = 0; period = 2147483646_8; gname = 'MINSTD2'
  case (5); M = 4294967296_8; a = 69069_8; c = 1; period = 4294967296_8; gname = 'LCG69069'
  case default; stop 'bad gen'
  end select
  if (nsteps < 0) nsteps = period

  ! ---- cross-check the local LCG against rnglib for the course generators
  if (gen <= 3) then
     ng = G_RANDU; if (gen == 2) ng = G_RANDM; if (gen == 3) ng = G_MINSTD
     call rng_init(ng, x0)
     x = x0
     do i = 1, 1000000
        snorm = rng_next(ng)
        x = step(x)
        if (x /= rng_state_int(ng)) stop 'local LCG differs from rnglib'
     end do
     print '(a,a)', trim(gname), ': local LCG matches rnglib for 1e6 steps'
  end if

  ! ---- read the vectors for this generator
  open(10, file='vectors.txt', status='old')
  nv = 0
  do
     read(10, '(a)', iostat=ios) line
     if (ios /= 0) exit
     read(line, *) ng
     if (ng /= gen) cycle
     nv = nv + 1
     read(line, *) ng, tv(nv), mode(nv), gran(nv), hnum(nv), hden(nv), s(:, nv)
  end do
  close(10)

  ! plane-index ranges: r/M lies in (sum of negative s_i, sum of positive s_i)
  l1 = 0
  do iv = 1, nv
     jlo(iv) = sum(min(s(1:tv(iv), iv), 0_8)) - 1
     jhi(iv) = sum(max(s(1:tv(iv), iv), 0_8)) + 1
     l1 = max(l1, jhi(iv) - jlo(iv) + 1)
  end do
  ! for mode 2 the index is r / gran, range (jhi-jlo)*M/gran
  do iv = 1, nv
     if (mode(iv) == 2) l1 = max(l1, (jhi(iv) - jlo(iv)) * (M / gran(iv)) + 2)
  end do
  allocate(cnt(0:l1, MAXRES, nv))
  cnt = 0; nres = 0; overflow = .false.

  ! ---- fill the window with the first 7 states, then run nsteps tuples
  x = x0
  win = 0
  do i = 2, 8
     win(i) = x
     x = step(x)
  end do
  ! now win(2..8) = X_0..X_6 and x = X_7
  x7 = x
  halfret = .false.
  call cpu_time(t0)
  do n = 1, nsteps
     win(1:7) = win(2:8)
     win(8) = x
     x = step(x)
     if (n == nsteps / 2) halfret = (x == x7)
     do iv = 1, nv
        t = tv(iv)
        r = 0
        do i = 1, t
           r = r + s(i, iv) * win(8 - t + i)
        end do
        if (mode(iv) == 1) then
           rr = modulo(r, M)
           j = (r - rr) / M
           do k = 1, int(nres(iv))
              if (res(k, iv) == rr) exit
           end do
           if (k > nres(iv)) then
              if (nres(iv) < MAXRES) then
                 nres(iv) = nres(iv) + 1; res(nres(iv), iv) = rr; k = int(nres(iv))
              else
                 overflow(iv) = .true.; cycle
              end if
           end if
           cnt(j - jlo(iv), k, iv) = cnt(j - jlo(iv), k, iv) + 1
        else
           if (modulo(r, gran(iv)) /= 0) overflow(iv) = .true.
           ir = (r - jlo(iv) * M) / gran(iv)
           cnt(ir, 1, iv) = cnt(ir, 1, iv) + 1
        end if
     end do
  end do
  call cpu_time(t1)
  print '(a,a,i12,a,f8.1,a)', trim(gname), ' tuples checked: ', nsteps, '  cpu ', t1 - t0, ' s'
  if (nsteps == period) print '(a,l2,a,l2)', '  state after P steps equals start:', (x == x7), &
       '   state after P/2 steps equals start:', halfret

  ! ---- report
  open(11, file='verify_'//trim(gname)//'.txt', status='replace')
  write(11, '(a)') '# t mode h  s(1:t) | nres residues | planes jmin jmax | emp_maxgap pred_maxgap ratio'
  do iv = 1, nv
     t = tv(iv)
     snorm = sqrt(dble(sum(s(1:t, iv)**2)))
     ! sorted plane positions
     np = 0
     allocate(pos(count(cnt(:, :, iv) > 0)))
     if (mode(iv) == 1) then
        do j = 0, l1
           do kk = 1, MAXRES   ! residues ascending within the same j
              do k = 1, int(nres(iv))
                 if (rank_of(res(1:nres(iv), iv), k) /= kk) cycle
                 if (cnt(j, k, iv) > 0) then
                    np = np + 1; pos(np) = (j + jlo(iv)) * M + res(k, iv)
                 end if
              end do
           end do
        end do
     else
        do j = 0, l1
           if (cnt(j, 1, iv) > 0) then
              np = np + 1; pos(np) = jlo(iv) * M + j * gran(iv)
           end if
        end do
     end if
     gapmax = 0
     do k = 2, np
        gapmax = max(gapmax, dble(pos(k) - pos(k-1)))
     end do
     gapu = gapmax / (dble(M) * snorm)
     write(11, '(i2,i3,1x,i0,"/",i0,1x,8i7,a,i3,8(1x,i0))', advance='no') t, mode(iv), hnum(iv), hden(iv), &
          s(:, iv), ' | ', nres(iv), res(1:nres(iv), iv)
     write(11, '(a,i8,2i8,a,es14.6,es14.6,f10.6,a,l2)') ' | ', np, (pos(1) - modulo(pos(1), M)) / M, &
          (pos(np) - modulo(pos(np), M)) / M, ' | ', gapu, dble(hnum(iv)) / dble(hden(iv)) / snorm, &
          gapu / (dble(hnum(iv)) / dble(hden(iv)) / snorm), ' overflow', overflow(iv)
     ! per-plane populations for small plane counts
     if (np <= 64 .and. mode(iv) == 1) then
        write(11, '(a)', advance='no') '   populations:'
        do j = 0, l1
           do k = 1, int(nres(iv))
              if (cnt(j, k, iv) > 0) write(11, '(1x,i0,":",i0)', advance='no') j + jlo(iv), cnt(j, k, iv)
           end do
        end do
        write(11, *)
     end if
     deallocate(pos)
  end do
  close(11)
  print '(a)', 'wrote verify_'//trim(gname)//'.txt'

contains

  function step(xx) result(y)
    integer(8), intent(in) :: xx
    integer(8) :: y
    if (gen <= 2) then
       y = iand(a * xx, M - 1)
    else if (gen <= 4) then
       y = modulo(a * xx, M)
    else
       y = iand(a * xx + c, M - 1)
    end if
  end function step

  integer function rank_of(v, k)
    integer(8), intent(in) :: v(:)
    integer, intent(in) :: k
    integer :: q
    rank_of = 1
    do q = 1, size(v)
       if (v(q) < v(k)) rank_of = rank_of + 1
    end do
  end function rank_of

end program lat_verify
