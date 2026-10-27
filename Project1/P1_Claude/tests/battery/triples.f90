! =====================================================================
!  triples.f90 : all eight Legendre triple products at fixed lags
!     T_pqr(j,k) = < P_p(x_{n-k}) P_q(x_{n-j}) P_r(x_n) >,  p,q,r in {1,2}
!  (p on the EARLIEST value, r on the newest).  Null: mean 0, sd 1/sqrt(N).
!  Also accumulates the second moment of each product so the actual
!  standard error can be quoted under the alternative.
!  usage: triples gen seed log2N j k outfile
! =====================================================================
program triples
  use rnglib
  implicit none
  integer, parameter :: RB = 1024, MSK = RB - 1
  integer :: gen, lg, j, k, p, q, r, c
  integer(8) :: seed, NN, it, i0, ij, ik
  real(8) :: PV(2, 0:RB-1), x, t, s1(8), s2(8)
  real(8), parameter :: R3 = 1.7320508075688772d0, R5 = 2.2360679774997897d0
  character(len=256) :: arg, outfile
  call get_command_argument(1, arg); read(arg,*) gen
  call get_command_argument(2, arg); read(arg,*) seed
  call get_command_argument(3, arg); read(arg,*) lg
  call get_command_argument(4, arg); read(arg,*) j
  call get_command_argument(5, arg); read(arg,*) k
  call get_command_argument(6, outfile)
  NN = 2_8**lg
  call rng_init(gen, seed)
  s1 = 0; s2 = 0
  do it = 0, NN - 1
     x = rng_next(gen)
     i0 = iand(it, int(MSK,8))
     PV(1, i0) = R3 * (2.0d0*x - 1.0d0)
     PV(2, i0) = R5 * (6.0d0*x*x - 6.0d0*x + 1.0d0)
     if (it >= k) then
        ij = iand(it - j, int(MSK,8)); ik = iand(it - k, int(MSK,8))
        c = 0
        do p = 1, 2
           do q = 1, 2
              do r = 1, 2
                 c = c + 1
                 t = PV(p, ik) * PV(q, ij) * PV(r, i0)
                 s1(c) = s1(c) + t
                 s2(c) = s2(c) + t*t
              end do
           end do
        end do
     end if
  end do
  open(11, file=trim(outfile), status='replace')
  write(11,'(a,a,a,i12,a,i14,a,2i5)') '# ', gen_name(gen), ' seed ', seed, ' N ', NN, ' lags j k ', j, k
  write(11,'(a)') '# p q r  mean  mean_sq   (p on x_{n-k}, q on x_{n-j}, r on x_n)'
  c = 0
  do p = 1, 2
     do q = 1, 2
        do r = 1, 2
           c = c + 1
           write(11,'(3i3,2es24.15)') p, q, r, s1(c)/dble(NN-k), s2(c)/dble(NN-k)
        end do
     end do
  end do
  close(11)
end program
