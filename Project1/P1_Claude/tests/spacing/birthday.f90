! =====================================================================
!  birthday.f90 : Marsaglia's birthday-spacings test in t = 1, 2, 3
!  dimensions (Marsaglia 1985; Knuth TAOCP 3.3.2 J; L'Ecuyer & Simard
!  2001), a "spacings" test that also sees multidimensional structure.
!
!  One replicate: n points in [0,1)^t, each built from t consecutive
!  outputs (non-overlapping); cell index c = sum_l floor(d u_l) d^(l-1)
!  with d = 2^b, so k = d^t "days". Sort the n birthdays, form the n-1
!  spacings between consecutive sorted birthdays (plus the first
!  birthday itself, as in Knuth), sort the spacings, and count
!      Y = number of spacing values equal to the previous one.
!  Under H0, Y is approximately Poisson with lambda = n^3/(4k).
!
!  usage: birthday gen seed single planfile outfile
!  plan lines: t b n R     output lines: t b n lambda r Y
! =====================================================================
program birthday
  use rnglib
  implicit none
  integer :: gen, single, ios, t, b, r, rr, l, i, nb
  integer(8) :: seed, n, y, d, c
  integer(8), allocatable :: bd(:), sp(:)
  real(8) :: lam
  character(len=256) :: arg, planfile, outfile
  integer :: pt(100), pb(100), pr(100)
  integer(8) :: pn(100)

  call get_command_argument(1, arg); read(arg, *) gen
  call get_command_argument(2, arg); read(arg, *) seed
  call get_command_argument(3, arg); read(arg, *) single
  call get_command_argument(4, planfile)
  call get_command_argument(5, outfile)
  rng_single = (single /= 0)
  open(10, file=planfile, status='old')
  nb = 0
  do
     read(10, *, iostat=ios) t, b, n, r
     if (ios /= 0) exit
     nb = nb + 1
     pt(nb) = t; pb(nb) = b; pn(nb) = n; pr(nb) = r
  end do
  close(10)

  call rng_init(gen, seed)
  open(20, file=outfile, status='replace')
  write(20, '(a,a,a,i0,a,i0)') '# ', trim(gen_name(gen)), ' seed ', seed, ' single ', single
  write(20, '(a)') '# t b n lambda r Y'
  do i = 1, nb
     t = pt(i); b = pb(i); n = pn(i)
     d = 2_8**b
     lam = dble(n)**3 / (4.0d0 * 2.0d0**(t * b))
     allocate(bd(n), sp(n))
     do rr = 1, pr(i)
        do c = 1, n
           bd(c) = 0
           do l = 1, t
              bd(c) = bd(c) * d + min(int(rng_next(gen) * dble(d), 8), d - 1)
           end do
        end do
        call isort(bd, n)
        sp(1) = bd(1)
        do c = 2, n
           sp(c) = bd(c) - bd(c - 1)
        end do
        call isort(sp, n)
        y = 0
        do c = 2, n
           if (sp(c) == sp(c - 1)) y = y + 1
        end do
        write(20, '(i0,1x,i0,1x,i0,1x,es14.6,1x,i0,1x,i0)') t, b, n, lam, rr, y
     end do
     deallocate(bd, sp)
  end do
  close(20)

contains

  subroutine isort(a, n)      ! heap sort of integer(8) a(1:n)
    integer(8), intent(inout) :: a(:)
    integer(8), intent(in) :: n
    integer(8) :: i, last, tmp
    do i = n / 2, 1, -1
       call isift(a, i, n)
    end do
    do last = n, 2, -1
       tmp = a(1); a(1) = a(last); a(last) = tmp
       call isift(a, 1_8, last - 1)
    end do
  end subroutine isort

  subroutine isift(a, start, m)
    integer(8), intent(inout) :: a(:)
    integer(8), intent(in) :: start, m
    integer(8) :: k, cc, tmp
    k = start
    tmp = a(k)
    do
       cc = 2 * k
       if (cc > m) exit
       if (cc < m) then
          if (a(cc + 1) > a(cc)) cc = cc + 1
       end if
       if (tmp >= a(cc)) exit
       a(k) = a(cc)
       k = cc
    end do
    a(k) = tmp
  end subroutine isift

end program birthday
