! =====================================================================
! lat_samples.f90 : compact samples for the lattice figures
!  (1) pairs (x_k, x_{k+1}) with x_k < 3/6539 among 1e8 overlapping pairs,
!      for RANDM, RANDU, MT19937               -> pairs_origin_<gen>.txt
!  (2) edge-on projections of non-overlapping t-tuples onto the shortest
!      dual vector s:  u = s.x/|s|,  v = w.x with w a unit vector orthogonal
!      to s (the diagonal (1,...,1) with its s-component removed)
!                                               -> edgeon_<tag>.txt
!  (3) RANDM triplets in a thin slab of the projection onto (25,115,678)
!                                               -> randm3_slab_<gen>.txt
! =====================================================================
program lat_samples
  use rnglib
  implicit none
  integer, parameter :: NE = 7
  integer(8) :: k, nsel, st4, st5
  real(8) :: x, xp, lim, y(8), s(8), w(8), u, v, sn
  integer :: ig, ids(3), unit, i, t, ie, gens(NE), tt(NE), nkeep
  character(len=12) :: nm3(3), tag(NE)
  integer :: sv(8, NE)

  ! ---------------- (1) pairs near the origin
  ids = [G_RANDM, G_RANDU, G_MT]
  nm3 = [character(len=12) :: 'RANDM', 'RANDU', 'MT19937']
  lim = 3.0d0 / 6539.0d0
  do ig = 1, 3
     call rng_init(ids(ig), 777_8)
     open(newunit=unit, file='pairs_origin_'//trim(nm3(ig))//'.txt', status='replace')
     xp = rng_next(ids(ig)); nsel = 0
     do k = 1, 100000000_8
        x = rng_next(ids(ig))
        if (xp < lim) then
           write(unit, '(es16.9,1x,f12.9)') xp, x
           nsel = nsel + 1
        end if
        xp = x
     end do
     close(unit)
     print '(a,a,i8)', trim(nm3(ig)), ' pairs with x_k < 3/6539 out of 1e8: ', nsel
  end do

  ! ---------------- (2) edge-on gallery
  ! gen codes: 1 RANDU 2 RANDM 3 MINSTD 4 MINSTD2 5 LCG69069 6 MT19937
  gens = [1, 2, 3, 4, 5, 6, 6]
  tt   = [3, 6, 8, 8, 6, 3, 6]
  tag  = [character(len=12) :: 'RANDU_t3', 'RANDM_t6', 'MINSTD_t8', 'MINSTD2_t8', 'LCG69069_t6', &
          'MT_t3', 'MT_t6']
  sv = 0
  sv(1:3, 1) = [9, -6, 1]
  sv(1:6, 2) = [3, 6, -8, -14, 9, -4]
  sv(1:8, 3) = [5, 4, 8, 1, 7, -2, 0, 1]
  sv(1:8, 4) = [3, -4, 4, 4, 2, -1, 2, 4]
  sv(1:6, 5) = [3, -11, 4, -8, 4, -4]
  sv(1:3, 6) = [9, -6, 1]
  sv(1:6, 7) = [3, 6, -8, -14, 9, -4]
  do ie = 1, NE
     t = tt(ie)
     s = 0; s(1:t) = dble(sv(1:t, ie)); sn = sqrt(sum(s * s)); s = s / sn
     w = 0; w(1:t) = 1.0d0
     w(1:t) = w(1:t) - sum(w(1:t) * s(1:t)) * s(1:t)
     w = w / sqrt(sum(w * w))
     call ginit(gens(ie))
     open(newunit=unit, file='edgeon_'//trim(tag(ie))//'.txt', status='replace')
     write(unit, '(a,8i6)') '# u = s.x/|s|, v = w.x ; s =', sv(1:t, ie)
     do k = 1, 20000
        do i = 1, t
           y(i) = gdraw(gens(ie))
        end do
        u = sum(s(1:t) * y(1:t)); v = sum(w(1:t) * y(1:t))
        write(unit, '(2f12.8)') u, v
     end do
     close(unit)
  end do

  ! ---------------- (3) RANDM 3D slab (8 planes wide)
  s = 0; s(1:3) = [25d0, 115d0, 678d0]; sn = sqrt(sum(s * s)); s = s / sn
  w = 0; w(1:3) = [0d0, 678d0, -115d0]; w = w / sqrt(sum(w * w))   ! orthogonal to s
  do ig = 1, 2
     if (ig == 1) call rng_init(G_RANDM, 777_8)
     if (ig == 2) call rng_init(G_MT, 777_8)
     open(newunit=unit, file='randm3_slab_'//trim(merge('RANDM  ', 'MT19937', ig == 1))//'.txt', status='replace')
     nkeep = 0
     do k = 1, 3000000
        y(1) = rng_next(merge(G_RANDM, G_MT, ig == 1))
        y(2) = rng_next(merge(G_RANDM, G_MT, ig == 1))
        y(3) = rng_next(merge(G_RANDM, G_MT, ig == 1))
        u = sum(s(1:3) * y(1:3))
        if (u >= 0.60d0 .and. u < 0.60d0 + 8.0d0 / sn) then
           v = sum(w(1:3) * y(1:3))
           write(unit, '(2f13.9)') u, v
           nkeep = nkeep + 1
        end if
     end do
     close(unit)
     print '(a,i2,a,i8)', 'slab sample gen', ig, ' kept ', nkeep
  end do

contains

  subroutine ginit(g)
    integer, intent(in) :: g
    select case (g)
    case (1); call rng_init(G_RANDU, 777_8)
    case (2); call rng_init(G_RANDM, 777_8)
    case (3); call rng_init(G_MINSTD, 777_8)
    case (4); st4 = 777_8
    case (5); st5 = 777_8
    case (6); call rng_init(G_MT, 777_8)
    end select
  end subroutine ginit

  function gdraw(g) result(r)
    integer, intent(in) :: g
    real(8) :: r
    select case (g)
    case (1); r = rng_next(G_RANDU)
    case (2); r = rng_next(G_RANDM)
    case (3); r = rng_next(G_MINSTD)
    case (4); st4 = modulo(48271_8 * st4, 2147483647_8); r = dble(st4) / 2147483647d0
    case (5); st5 = iand(69069_8 * st5 + 1_8, 4294967295_8); r = dble(st5) / 4294967296d0
    case default; r = rng_next(G_MT)
    end select
  end function gdraw

end program lat_samples
