! dump window values (and global draw index) of replicate rep of block (N, j)
program dumpwin
  use rnglib
  implicit none
  integer :: gen, j, rep, r
  integer(8) :: seed, nn, i, n
  real(8) :: u, lo
  character(len=256) :: arg, outf
  real(8), allocatable :: v(:)
  integer(8), allocatable :: idx(:)
  call get_command_argument(1, arg); read(arg,*) gen
  call get_command_argument(2, arg); read(arg,*) seed
  call get_command_argument(3, arg); read(arg,*) nn
  call get_command_argument(4, arg); read(arg,*) j
  call get_command_argument(5, arg); read(arg,*) rep
  call get_command_argument(6, outf)
  lo = 1.0d0 - 2.0d0**(-j)
  allocate(v(nn/2**j*2+1000), idx(nn/2**j*2+1000))
  call rng_init(gen, seed)
  do r = 1, rep
    n = 0
    do i = 1, nn
      u = rng_next(gen)
      if (u >= lo) then
        n = n + 1
        if (r == rep) then
          v(n) = u; idx(n) = (r-1)*nn + i
        end if
      end if
    end do
  end do
  open(30, file=outf, access='stream', form='unformatted', status='replace')
  write(30) n
  write(30) v(1:n)
  write(30) idx(1:n)
  close(30)
end program dumpwin
