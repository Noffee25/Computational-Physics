! =====================================================================
!  spacing_scan.f90 : maximum spacing, minimum spacing and duplicate
!  counts for one generator configuration, over a plan of (N, R, j).
!
!  For every block (N, R, j) of the plan the generator is re-seeded and
!  R consecutive replicates of N draws are taken. Only draws in the top
!  window [1-a, 1), a = 2^-j, are kept (j = 0 keeps everything). The
!  kept points are rescaled exactly to v = (u - (1-a)) * 2^j in [0,1)
!  (Sterbenz: the subtraction is exact for a <= 1/2) and sorted with an
!  O(n) bucket sort. With the endpoints 0 and 1 the n points define
!  n+1 spacings D_0..D_n. Per replicate we write
!      n, M = max D_i, m = min D_i, dup = n - #distinct values.
!
!  usage: spacing_scan gen seed single planfile outfile
! =====================================================================
program spacing_scan
  use rnglib
  implicit none
  integer :: gen, single, ios, nplan, ib, r, j, flag
  integer(8) :: seed, nn, i, n, cap, capmax, dup
  integer(8), allocatable :: plan_n(:)
  integer, allocatable :: plan_r(:), plan_j(:)
  real(8), allocatable :: v(:), s(:)
  integer, allocatable :: cnt(:)
  real(8) :: u, a, lo, sc, bigm, smallm
  character(len=256) :: arg, planfile, outfile
  integer(8) :: t0, t1, rate

  call get_command_argument(1, arg); read(arg, *) gen
  call get_command_argument(2, arg); read(arg, *) seed
  call get_command_argument(3, arg); read(arg, *) single
  call get_command_argument(4, planfile)
  call get_command_argument(5, outfile)
  rng_single = (single /= 0)

  ! ---- read plan
  open(10, file=planfile, status='old')
  nplan = 0
  do
     read(10, *, iostat=ios) nn, r, j
     if (ios /= 0) exit
     nplan = nplan + 1
  end do
  allocate(plan_n(nplan), plan_r(nplan), plan_j(nplan))
  rewind(10)
  do ib = 1, nplan
     read(10, *) plan_n(ib), plan_r(ib), plan_j(ib)
  end do
  close(10)

  capmax = 16000000_8
  open(20, file=outfile, status='replace')
  write(20, '(a)') '# gen seed single : ' // trim(gen_name(gen))
  write(20, '(a,i0,a,i0,a,i0)') '# ', gen, ' ', seed, ' ', single
  write(20, '(a)') '# N R_index j n M m dup flag'
  call system_clock(t0, rate)

  do ib = 1, nplan
     nn = plan_n(ib); j = plan_j(ib)
     a  = 2.0d0**(-j)
     lo = 1.0d0 - a
     sc = 2.0d0**j
     cap = min(capmax, int(dble(nn) * a * 1.02d0 + 6.0d0 * sqrt(dble(nn) * a) + 100.0d0, 8))
     if (j == 0) cap = min(capmax, nn)
     allocate(v(cap), s(cap), cnt(0:cap))
     call rng_init(gen, seed)
     do r = 1, plan_r(ib)
        n = 0
        flag = 0
        if (j == 0) then
           do i = 1, nn
              u = rng_next(gen)
              n = n + 1
              if (n <= cap) v(n) = u
           end do
        else
           do i = 1, nn
              u = rng_next(gen)
              if (u >= lo) then
                 n = n + 1
                 if (n <= cap) v(n) = (u - lo) * sc
              end if
           end do
        end if
        if (n > cap) then
           flag = 2            ! window overflow: far more points than expected
           bigm = -1.0d0; smallm = -1.0d0; dup = -1
        else
           call bucket_sort(v, s, cnt, n)
           call spacings(s, n, bigm, smallm, dup)
        end if
        write(20, '(i0,1x,i0,1x,i0,1x,i0,1x,es24.16e3,1x,es24.16e3,1x,i0,1x,i0)') &
             nn, r, j, n, bigm, smallm, dup, flag
     end do
     deallocate(v, s, cnt)
     call system_clock(t1)
     write(*, '(a,a,a,i0,a,i0,a,f8.1,a)') trim(gen_name(gen)), ' ', 'N=', nn, ' R=', plan_r(ib), &
          ' t=', dble(t1 - t0) / dble(rate), ' s'
     flush(20)
  end do
  close(20)

contains

  ! O(n) bucket sort of v(1:n) in [0,1) into s(1:n)
  subroutine bucket_sort(v, s, cnt, n)
    real(8), intent(in)    :: v(:)
    real(8), intent(inout) :: s(:)
    integer, intent(inout) :: cnt(0:)
    integer(8), intent(in) :: n
    integer :: nb, b, k, lo_, hi_
    integer(8) :: i
    if (n <= 0) return
    nb = int(n)
    cnt(0:nb) = 0
    do i = 1, n
       b = min(int(v(i) * nb), nb - 1)
       cnt(b + 1) = cnt(b + 1) + 1
    end do
    do k = 1, nb
       cnt(k) = cnt(k) + cnt(k - 1)
    end do
    ! cnt(b) is now the start offset (0-based) of bucket b
    do i = 1, n
       b = min(int(v(i) * nb), nb - 1)
       cnt(b) = cnt(b) + 1
       s(cnt(b)) = v(i)
    end do
    ! after placement cnt(b) = end index (1-based) of bucket b; bucket b
    ! occupies (cnt(b-1)+1 .. cnt(b)) with cnt(-1) = 0
    lo_ = 1
    do b = 0, nb - 1
       hi_ = cnt(b)
       if (hi_ - lo_ >= 1) then
          if (hi_ - lo_ < 32) then
             call insertion_sort(s, lo_, hi_)
          else
             call heap_sort(s, lo_, hi_)
          end if
       end if
       lo_ = hi_ + 1
    end do
  end subroutine bucket_sort

  subroutine insertion_sort(s, lo_, hi_)
    real(8), intent(inout) :: s(:)
    integer, intent(in) :: lo_, hi_
    integer :: i, k
    real(8) :: t
    do i = lo_ + 1, hi_
       t = s(i)
       k = i - 1
       do while (k >= lo_)
          if (s(k) <= t) exit
          s(k + 1) = s(k)
          k = k - 1
       end do
       s(k + 1) = t
    end do
  end subroutine insertion_sort

  subroutine heap_sort(s, lo_, hi_)
    real(8), intent(inout) :: s(:)
    integer, intent(in) :: lo_, hi_
    integer :: n, i
    real(8) :: t
    n = hi_ - lo_ + 1
    do i = n / 2, 1, -1
       call sift(s, lo_, i, n)
    end do
    do i = n, 2, -1
       t = s(lo_); s(lo_) = s(lo_ + i - 1); s(lo_ + i - 1) = t
       call sift(s, lo_, 1, i - 1)
    end do
  end subroutine heap_sort

  subroutine sift(s, lo_, start, last)
    real(8), intent(inout) :: s(:)
    integer, intent(in) :: lo_, start, last
    integer :: k, c
    real(8) :: t
    k = start
    t = s(lo_ + k - 1)
    do
       c = 2 * k
       if (c > last) exit
       if (c < last) then
          if (s(lo_ + c) > s(lo_ + c - 1)) c = c + 1
       end if
       if (t >= s(lo_ + c - 1)) exit
       s(lo_ + k - 1) = s(lo_ + c - 1)
       k = c
    end do
    s(lo_ + k - 1) = t
  end subroutine sift

  subroutine spacings(s, n, bigm, smallm, dup)
    real(8), intent(in) :: s(:)
    integer(8), intent(in) :: n
    real(8), intent(out) :: bigm, smallm
    integer(8), intent(out) :: dup
    integer(8) :: i
    real(8) :: d
    if (n <= 0) then
       bigm = 1.0d0; smallm = 1.0d0; dup = 0; return
    end if
    bigm = s(1); smallm = s(1); dup = 0
    do i = 2, n
       d = s(i) - s(i - 1)
       if (d > bigm) bigm = d
       if (d < smallm) smallm = d
       if (d == 0.0d0) dup = dup + 1
    end do
    d = 1.0d0 - s(n)
    if (d > bigm) bigm = d
    if (d < smallm) smallm = d
  end subroutine spacings

end program spacing_scan
