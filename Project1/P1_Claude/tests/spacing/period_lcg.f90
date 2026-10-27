! =====================================================================
!  period_lcg.f90 : exact period of RANDU, RANDM (x -> a x mod 2^31)
!  and MINSTD (x -> 16807 x mod 2^31-1) by Brent's cycle detection on
!  the full integer state, for a list of seeds.
!
!  The inner maps are re-implemented with fast integer arithmetic; the
!  first 10^6 states are checked against rnglib (rng_state_int) so the
!  measured periods are those of the shared library generators.
!
!  output lines: gen seed trailing_zeros lambda mu seconds
! =====================================================================
program period_lcg
  use rnglib
  implicit none
  integer(8), parameter :: MASK31 = 2147483647_8
  integer(8) :: seeds(64), lam, mu, x, y
  integer :: nseed, i, k, g, gens(3)
  integer(8) :: t0, t1, rate
  integer(8), parameter :: odd(7) = [3_8, 77_8, 777_8, 7777_8, 12345_8, 77777_8, 314159_8]

  ! ---- verify fast maps against rnglib
  gens = [G_RANDU, G_RANDM, G_MINSTD]
  do g = 1, 3
     call rng_init(gens(g), 777_8)
     x = 777_8
     do i = 1, 1000000
        y = rng_next(gens(g))
        x = fmap(gens(g), x)
        if (x /= rng_state_int(gens(g))) then
           print *, 'MISMATCH ', gen_name(gens(g)), i, x, rng_state_int(gens(g))
           stop 1
        end if
     end do
     print '(a,a)', '# fast map verified against rnglib for ', trim(gen_name(gens(g)))
  end do

  ! ---- seed list: odd seeds with 1..6 digits, 2^k*777 and 2^k (k=1..24), 0
  nseed = 0
  do i = 1, 7
     nseed = nseed + 1; seeds(nseed) = odd(i)
  end do
  do k = 1, 24
     nseed = nseed + 1; seeds(nseed) = iand(777_8 * 2_8**k, MASK31)
  end do
  do k = 1, 24
     nseed = nseed + 1; seeds(nseed) = 2_8**k
  end do
  nseed = nseed + 1; seeds(nseed) = 0_8

  call system_clock(t0, rate)
  open(20, file='period_lcg.txt', status='replace')
  write(20, '(a)') '# gen seed tz lambda mu seconds'
  do g = 1, 2
     do i = 1, nseed
        call system_clock(t0)
        call brent(gens(g), seeds(i), lam, mu)
        call system_clock(t1)
        write(20, '(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,f8.2)') trim(gen_name(gens(g))), seeds(i), &
             tz(seeds(i)), lam, mu, dble(t1 - t0) / dble(rate)
        flush(20)
     end do
  end do
  ! MINSTD: seeds 1 and 777 (16807 is a primitive root mod 2^31-1)
  do i = 1, 2
     call system_clock(t0)
     if (i == 1) then
        call brent(G_MINSTD, 1_8, lam, mu)
     else
        call brent(G_MINSTD, 777_8, lam, mu)
     end if
     call system_clock(t1)
     write(20, '(a,1x,i0,1x,i0,1x,i0,1x,i0,1x,f8.2)') 'MINSTD', merge(1_8, 777_8, i == 1), 0, lam, mu, &
          dble(t1 - t0) / dble(rate)
     flush(20)
  end do
  close(20)

contains

  pure function fmap(gid, x) result(y)
    integer, intent(in) :: gid
    integer(8), intent(in) :: x
    integer(8) :: y
    select case (gid)
    case (G_RANDU);  y = iand(65539_8 * x, MASK31)
    case (G_RANDM);  y = iand(6539_8 * x, MASK31)
    case default      ! MINSTD, Mersenne reduction mod 2^31-1
       y = 16807_8 * x
       y = iand(y, MASK31) + ishft(y, -31)
       if (y >= MASK31) y = y - MASK31
    end select
  end function fmap

  pure function tz(s) result(k)
    integer(8), intent(in) :: s
    integer :: k
    if (s == 0) then
       k = 31
    else
       k = trailz(s)
    end if
  end function tz

  ! Brent's algorithm: lambda = cycle length, mu = transient length
  subroutine brent(gid, x0, lam, mu)
    integer, intent(in) :: gid
    integer(8), intent(in) :: x0
    integer(8), intent(out) :: lam, mu
    integer(8) :: power, tort, hare, i
    power = 1; lam = 1
    tort = x0; hare = fmap(gid, x0)
    do while (tort /= hare)
       if (power == lam) then
          tort = hare; power = 2 * power; lam = 0
       end if
       hare = fmap(gid, hare); lam = lam + 1
    end do
    tort = x0; hare = x0
    do i = 1, lam
       hare = fmap(gid, hare)
    end do
    mu = 0
    do while (tort /= hare)
       tort = fmap(gid, tort); hare = fmap(gid, hare); mu = mu + 1
    end do
  end subroutine brent

end program period_lcg
