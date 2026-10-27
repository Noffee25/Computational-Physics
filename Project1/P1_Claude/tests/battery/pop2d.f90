! pop2d.f90 : exact 2D cell populations of non-overlapping pairs over one
! full period, for RANDM (2^28 pairs) or MINSTD (2^30-1 pairs), k=1000.
! Prints relvar = C * sum (p_i - 1/C)^2 and the period of the pair sequence.
program pop2d
  use rnglib
  implicit none
  integer, parameter :: k = 1000
  integer, allocatable :: h(:,:)
  integer :: gen
  integer(8) :: s0, s, t, np
  real(8) :: x, y, C, rv, p
  character(len=16) :: arg
  integer :: i, j
  call get_command_argument(1, arg); read(arg,*) gen
  allocate(h(0:k-1,0:k-1)); h = 0
  call rng_init(gen, 777_8)
  s0 = rng_state_int(gen)
  np = 0
  do
     x = rng_next(gen); y = rng_next(gen); np = np + 1
     h(min(int(x*k),k-1), min(int(y*k),k-1)) = h(min(int(x*k),k-1), min(int(y*k),k-1)) + 1
     if (rng_state_int(gen) == s0) exit
  end do
  C = dble(k)**2
  rv = 0
  do j = 0, k-1
     do i = 0, k-1
        p = dble(h(i,j)) / dble(np)
        rv = rv + (p - 1d0/C)**2
     end do
  end do
  rv = rv * C
  print '(a,a,i12,a,es14.6,a,es14.6)', gen_name(gen), ' pairs in period ', np, ' relvar ', rv, &
       '  N_t*relvar for 1e7 pairs ', 1d7*rv
end program
