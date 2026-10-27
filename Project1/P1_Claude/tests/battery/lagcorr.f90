! =====================================================================
!  lagcorr.f90 : lag-k Legendre correlations of a generator (GEMM form)
!
!  P1(x) = sqrt(3)(2x-1), P2(x) = sqrt(5)(6x^2-6x+1) are orthonormal
!  under U(0,1).  For lag k = 1..L we accumulate
!     S_pq(k) = sum_n P_p(x_n) P_q(x_{n-k})      (p on the LATER value)
!  so that r_pq(k) = S_pq(k)/(N-k) has mean 0 and sd 1/sqrt(N-k) for
!  an ideal generator.  r_11 is the usual serial correlation.
!
!  Method: the stream is cut into sub-blocks of length L, stored as the
!  columns of E(2L,0:M) (rows 1..L: P1, rows L+1..2L: P2).  Then
!     C = sum_c E(:,c) E(:,c)^T,   A = sum_c E(:,c) E(:,c-1)^T
!  (two DGEMMs) contain every product with lag <= L:
!     S_pq(k) = sum_{l>k} C_pq(l,l-k) + sum_{l<=k} A_pq(l,l+L-k).
!
!  usage: lagcorr gen seed log2N L outfile
! =====================================================================
program lagcorr
  use rnglib
  implicit none
  integer, parameter :: M = 4096
  integer :: gen, L, lg, k, j, c, ib, nb, L2
  integer(8) :: seed, N, done, ndone
  real(8), allocatable :: E(:,:), Ct(:,:), At(:,:)
  real(8), allocatable :: S(:,:)
  real(8) :: x, sx, sxx, t0, t1
  real(8), parameter :: R3 = 1.7320508075688772d0, R5 = 2.2360679774997897d0
  character(len=256) :: arg, outfile
  call get_command_argument(1, arg); read(arg,*) gen
  call get_command_argument(2, arg); read(arg,*) seed
  call get_command_argument(3, arg); read(arg,*) lg
  call get_command_argument(4, arg); read(arg,*) L
  call get_command_argument(5, outfile)
  N = 2_8**lg
  L2 = 2*L
  allocate(E(L2, 0:M), Ct(L2, L2), At(L2, L2), S(L, 4))
  E = 0; Ct = 0; At = 0
  sx = 0; sxx = 0
  call rng_init(gen, seed)
  nb = int((N + int(L,8)*M - 1) / (int(L,8)*M))
  done = 0
  open(10, file=trim(outfile)//'.cp', status='replace')   ! checkpoints for lag 1,2,3
  call cpu_time(t0)
  do ib = 1, nb
     E(:, 0) = E(:, M)
     do c = 1, M
        do j = 1, L
           if (done < N) then
              x = rng_next(gen)
              sx = sx + x; sxx = sxx + x*x
              E(j, c)     = R3 * (2.0d0*x - 1.0d0)
              E(L + j, c) = R5 * (6.0d0*x*x - 6.0d0*x + 1.0d0)
              done = done + 1
           else
              E(j, c) = 0; E(L + j, c) = 0
           end if
        end do
     end do
     call dgemm('N', 'T', L2, L2, M, 1.0d0, E(1,1), L2, E(1,1), L2, 1.0d0, Ct, L2)
     call dgemm('N', 'T', L2, L2, M, 1.0d0, E(1,1), L2, E(1,0), L2, 1.0d0, At, L2)
     ndone = done
     call extract(1, 1, 3)
     write(10,'(i12,3es24.15)') ndone, S(1,1), S(2,1), S(3,1)
  end do
  call cpu_time(t1)
  close(10)
  call extract(1, 1, L)
  call extract(1, 2, L)
  call extract(2, 1, L)
  call extract(2, 2, L)
  open(11, file=trim(outfile), status='replace')
  write(11,'(a,i3,a,a)') '# gen ', gen, ' ', gen_name(gen)
  write(11,'(a,i12)') '# seed ', seed
  write(11,'(a,i14)') '# N ', N
  write(11,'(a,2es24.15)') '# sumx sumxx ', sx, sxx
  write(11,'(a,f10.2)') '# cpu ', t1 - t0
  write(11,'(a)') '# k S11 S12 S21 S22   (S_pq = sum P_p(x_n) P_q(x_{n-k}); p on later)'
  do k = 1, L
     write(11,'(i5,4es24.15)') k, S(k,1), S(k,2), S(k,3), S(k,4)
  end do
  close(11)
  print '(a8,a,f8.2,a)', gen_name(gen), ' done in ', t1-t0, ' s'
contains
  subroutine extract(p, q, kmax)
    integer, intent(in) :: p, q, kmax
    integer :: kk, ll, op, oq, idx
    op = (p-1)*L; oq = (q-1)*L
    idx = 2*(p-1) + q          ! 11->1, 12->2, 21->3, 22->4
    do kk = 1, kmax
       S(kk, idx) = 0
       do ll = kk+1, L
          S(kk, idx) = S(kk, idx) + Ct(op+ll, oq+ll-kk)
       end do
       do ll = 1, kk
          S(kk, idx) = S(kk, idx) + At(op+ll, oq+ll+L-kk)
       end do
    end do
  end subroutine
end program
