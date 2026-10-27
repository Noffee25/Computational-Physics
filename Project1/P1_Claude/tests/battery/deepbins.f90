! deepbins.f90 : reproduce the student's 1e7-bin, 1e9-draw histogram and
! record the moments of the bin counts plus the count distribution.
! usage: deepbins gen seed ndraws nbins outfile
program deepbins
  use rnglib
  implicit none
  integer :: gen, nbins, idx, c, cmax
  integer(8) :: seed, nd, i, ncnt(0:400)
  integer, allocatable :: h(:)
  real(8) :: x, m1, m2
  character(len=256) :: arg, outfile
  call get_command_argument(1, arg); read(arg,*) gen
  call get_command_argument(2, arg); read(arg,*) seed
  call get_command_argument(3, arg); read(arg,*) nd
  call get_command_argument(4, arg); read(arg,*) nbins
  call get_command_argument(5, outfile)
  allocate(h(nbins)); h = 0
  call rng_init(gen, seed)
  do i = 1, nd
     x = rng_next(gen)
     idx = min(int(x * nbins) + 1, nbins)
     h(idx) = h(idx) + 1
  end do
  ncnt = 0; m1 = 0; m2 = 0; cmax = 0
  do idx = 1, nbins
     c = min(h(idx), 400)
     ncnt(c) = ncnt(c) + 1
     m1 = m1 + h(idx); m2 = m2 + dble(h(idx))**2
  end do
  m1 = m1 / nbins; m2 = m2 / nbins - m1*m1
  open(11, file=trim(outfile), status='replace')
  write(11,'(a,a,a,i12,a,i14,a,i10)') '# ', gen_name(gen), ' seed ', seed, ' draws ', nd, ' bins ', nbins
  write(11,'(a,2f14.6,2i8)') '# mean sd min max ', m1, sqrt(m2), minval(h), maxval(h)
  do c = 0, 400
     if (ncnt(c) > 0) write(11,'(i5,i12)') c, ncnt(c)
  end do
  close(11)
end program
