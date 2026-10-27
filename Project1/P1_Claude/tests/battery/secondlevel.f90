! =====================================================================
!  secondlevel.f90 : histogram chi-square on consecutive blocks
!  The stream (2^lgN values) is cut into blocks of 2^lgL values.  For each
!  block the K = 2^lgK bin histogram gives X^2_fine (K-1 df) and the
!  histogram coarsened to 2^10 bins gives X^2_coarse (1023 df).  The
!  cumulative histogram (prefix of length b*2^lgL) is also tested.
!  These first-level statistics are turned into p-values and tested for
!  uniformity (Kolmogorov-Smirnov) in battery_analysis.py.
!  usage: secondlevel gen seed lgN lgK lgL outfile
! =====================================================================
program secondlevel
  use rnglib
  implicit none
  integer :: gen, lgN, lgK, lgL, nblk, ib, KC, grp
  integer(8) :: seed, K, L, i, idx, j
  integer, allocatable :: h(:), hc(:)
  integer(8) :: s2, s2c, c2, c2c, csum
  real(8) :: x, xk, lam, lamc, t0, t1
  character(len=256) :: arg, outfile
  call get_command_argument(1, arg); read(arg,*) gen
  call get_command_argument(2, arg); read(arg,*) seed
  call get_command_argument(3, arg); read(arg,*) lgN
  call get_command_argument(4, arg); read(arg,*) lgK
  call get_command_argument(5, arg); read(arg,*) lgL
  call get_command_argument(6, outfile)
  K = 2_8**lgK; L = 2_8**lgL; KC = 1024; grp = int(K / KC)
  nblk = 2**(lgN - lgL)
  allocate(h(0:K-1), hc(0:K-1)); hc = 0
  xk = dble(K)
  call rng_init(gen, seed)
  open(11, file=trim(outfile), status='replace')
  write(11,'(a,a,a,i12,a,i3,a,i3,a,i3)') '# ', gen_name(gen), ' seed ', seed, ' lgN ', lgN, ' lgK ', lgK, ' lgL ', lgL
  write(11,'(a)') '# block  X2fine_block  X2coarse_block  X2fine_cum  X2coarse_cum'
  call cpu_time(t0)
  do ib = 1, nblk
     h = 0
     do i = 1, L
        x = rng_next(gen)
        idx = min(int(x * xk, 8), K - 1)
        h(idx) = h(idx) + 1
     end do
     hc = hc + h
     s2 = 0; c2 = 0; s2c = 0; c2c = 0
     do j = 0, K - 1
        s2 = s2 + int(h(j),8)**2
        c2 = c2 + int(hc(j),8)**2
     end do
     do j = 0, KC - 1
        csum = sum(int(h(j*grp:(j+1)*grp-1),8));  s2c = s2c + csum**2
        csum = sum(int(hc(j*grp:(j+1)*grp-1),8)); c2c = c2c + csum**2
     end do
     lam = dble(L) / dble(K); lamc = dble(L) / dble(KC)
     write(11,'(i6,4es24.15)') ib, dble(s2)/lam - dble(L), dble(s2c)/lamc - dble(L), &
          dble(c2)/(lam*ib) - dble(L)*ib, dble(c2c)/(lamc*ib) - dble(L)*ib
  end do
  call cpu_time(t1)
  close(11)
  print '(a8,a,f8.2)', gen_name(gen), ' secondlevel cpu ', t1-t0
end program
