! z-score of the hit-or-miss ball estimate after N = 1e7 triplets for many seeds
program lat_seeds
  use rnglib
  implicit none
  integer, parameter :: NS = 48
  integer(8), parameter :: N = 10000000_8
  integer(8) :: seeds(NS), nin, k, s
  integer :: is, ig, ids(2), unit
  real(8) :: v(3), p, q, z(NS, 2)
  ids = [G_RANDM, G_MT]
  q = 3.141592653589793d0 / 6d0
  do is = 1, 24
     seeds(is) = 2 * is - 1                  ! small odd seeds 1, 3, ..., 47
  end do
  s = 123456789_8
  do is = 25, NS                             ! large odd seeds from a simple LCG
     s = modulo(s * 1103515245_8 + 12345_8, 2147483648_8)
     seeds(is) = ior(s, 1_8)
  end do
  open(newunit=unit, file='ball_seeds.txt', status='replace')
  write(unit, '(a)') '# seed z_RANDM z_MT19937  (ball hit-or-miss, first 1e7 triplets)'
  do is = 1, NS
     do ig = 1, 2
        call rng_init(ids(ig), seeds(is))
        nin = 0
        do k = 1, N
           v(1) = 2d0 * rng_next(ids(ig)) - 1d0
           v(2) = 2d0 * rng_next(ids(ig)) - 1d0
           v(3) = 2d0 * rng_next(ids(ig)) - 1d0
           if (v(1)**2 + v(2)**2 + v(3)**2 <= 1d0) nin = nin + 1
        end do
        p = dble(nin) / dble(N)
        z(is, ig) = (p - q) / sqrt(q * (1 - q) / dble(N))
     end do
     write(unit, '(i12,2f9.3)') seeds(is), z(is, :)
  end do
  close(unit)
  print '(a,2f8.3)', 'mean z  RANDM, MT: ', sum(z(:, 1)) / NS, sum(z(:, 2)) / NS
  print '(a,2f8.3)', 'var  z  RANDM, MT: ', sum((z(:, 1) - sum(z(:, 1)) / NS)**2) / (NS - 1), &
       sum((z(:, 2) - sum(z(:, 2)) / NS)**2) / (NS - 1)
end program
