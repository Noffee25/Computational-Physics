! Hit-or-miss ball test in blocks: counts of triplets inside the unit ball
! (v = 2u - 1) for consecutive blocks of 1e4 non-overlapping triplets,
! over 178956970 triplets (= one full RANDU/RANDM period of draws).
program lat_blocks
  use rnglib
  implicit none
  integer, parameter :: NG = 7, B = 10000
  integer(8), parameter :: NT = 178956970_8
  integer :: ig, unit, nb, ib, i
  integer(8) :: k, st4, st5
  integer, allocatable :: cnt(:, :)
  real(8) :: v(3)
  character(len=8) :: nm(NG) = [character(len=8) :: 'RANDU', 'RANDM', 'MINSTD', 'MINSTD2', 'LCG69069', 'MT19937', 'RANNUM']
  nb = int(NT / B)
  allocate(cnt(nb, NG)); cnt = 0
  do ig = 1, NG
     select case (ig)
     case (1); call rng_init(G_RANDU, 777_8)
     case (2); call rng_init(G_RANDM, 777_8)
     case (3); call rng_init(G_MINSTD, 777_8)
     case (4); st4 = 777_8
     case (5); st5 = 777_8
     case (6); call rng_init(G_MT, 777_8)
     case (7); call rng_init(G_RANNUM, 777_8)
     end select
     do ib = 1, nb
        do k = 1, B
           do i = 1, 3
              v(i) = 2d0 * draw(ig) - 1d0
           end do
           if (v(1)**2 + v(2)**2 + v(3)**2 <= 1d0) cnt(ib, ig) = cnt(ib, ig) + 1
        end do
     end do
     print '(a,a,f12.8)', nm(ig), ' pi_hat over all blocks ', 6d0 * sum(dble(cnt(:, ig))) / (dble(nb) * B)
  end do
  open(newunit=unit, file='ball_blocks.txt', status='replace')
  write(unit, '(a,7a9)') '# hits per block of 1e4 triplets:', nm
  do ib = 1, nb
     write(unit, '(7i6)') cnt(ib, :)
  end do
  close(unit)
contains
  function draw(ig) result(y)
    integer, intent(in) :: ig
    real(8) :: y
    select case (ig)
    case (1); y = rng_next(G_RANDU)
    case (2); y = rng_next(G_RANDM)
    case (3); y = rng_next(G_MINSTD)
    case (4); st4 = modulo(48271_8 * st4, 2147483647_8); y = dble(st4) / 2147483647d0
    case (5); st5 = iand(69069_8 * st5 + 1_8, 4294967295_8); y = dble(st5) / 4294967296d0
    case (6); y = rng_next(G_MT)
    case default; y = rng_next(G_RANNUM)
    end select
  end function draw
end program
