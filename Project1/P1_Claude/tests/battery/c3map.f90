! =====================================================================
!  c3map.f90 : three-point correlation maps via DGEMM
!     T1(j,k) = < P1(x_n) P1(x_{n-j}) P1(x_{n-k}) >   (odd kernel; equals
!               12^{3/2} <(x_n-1/2)(x_{n-j}-1/2)(x_{n-k}-1/2)>)
!     T2(j,k) = < P2(x_n) P1(x_{n-j}) P1(x_{n-k}) >   (even total parity)
!  for 1 <= j,k <= K. Each entry has null mean 0 and sd 1/sqrt(N-max(j,k)).
!  For a block of B values, A(i,j) = P1(x_{i-j}) and W = [P1(x_i)A, P2(x_i)A],
!  so the block contribution is W^T A (one DGEMM, 4BK^2 flops).
!  usage: c3map gen seed log2N K outfile   (binary output: K, N, T(2K,K))
! =====================================================================
program c3map
  use rnglib
  implicit none
  integer, parameter :: B = 4096
  integer :: gen, K, lg, i, j, ib, nb
  integer(8) :: seed, N
  real(8), allocatable :: p1(:), p2(:), A(:,:), W(:,:), C(:,:)
  real(8) :: x, t0, t1
  real(8), parameter :: R3 = 1.7320508075688772d0, R5 = 2.2360679774997897d0
  character(len=256) :: arg, outfile
  call get_command_argument(1, arg); read(arg,*) gen
  call get_command_argument(2, arg); read(arg,*) seed
  call get_command_argument(3, arg); read(arg,*) lg
  call get_command_argument(4, arg); read(arg,*) K
  call get_command_argument(5, outfile)
  N = 2_8**lg
  nb = int(N / B)
  allocate(p1(1-K:B), p2(1:B), A(B,K), W(B,2*K), C(2*K,K))
  p1 = 0; C = 0
  call rng_init(gen, seed)
  call cpu_time(t0)
  do ib = 1, nb
     p1(1-K:0) = p1(B-K+1:B)
     do i = 1, B
        x = rng_next(gen)
        p1(i) = R3 * (2.0d0*x - 1.0d0)
        p2(i) = R5 * (6.0d0*x*x - 6.0d0*x + 1.0d0)
     end do
     do j = 1, K
        A(:, j) = p1(1-j:B-j)
        W(:, j) = p1(1:B) * A(:, j)
        W(:, K+j) = p2(1:B) * A(:, j)
     end do
     call dgemm('T', 'N', 2*K, K, B, 1.0d0, W, B, A, B, 1.0d0, C, 2*K)
  end do
  call cpu_time(t1)
  open(11, file=trim(outfile), access='stream', form='unformatted', status='replace')
  write(11) K, N, C
  close(11)
  print '(a8,a,f8.2,a)', gen_name(gen), ' c3map done in ', t1-t0, ' s'
end program
