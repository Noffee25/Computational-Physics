! =====================================================================
!  spacing_dups.f90 : cumulative duplicate count D(N) along one long
!  stream, for generators whose outputs live on a grid of at most 2^31
!  points, using an exact bitmap of the values already seen.
!
!  D(N) = N - (number of distinct values among the first N draws).
!  Keys: RANDU, RANDM, MINSTD use the integer state (equal to the output
!  times the modulus); RAND uses u*2^22 (exact integer); single
!  precision uses the IEEE bit pattern of real(u,4).
!
!  usage: spacing_dups gen seed single log2Nmax outfile
!  output lines: N D(N)
! =====================================================================
program spacing_dups
  use rnglib
  implicit none
  integer :: gen, single, l2max, i, ib
  integer(8) :: seed, n, nmax, key, w, d, nbits, nextcp
  integer(8), allocatable :: bm(:)
  real(8) :: u
  real(4) :: us
  character(len=256) :: arg, outfile
  integer :: icp

  call get_command_argument(1, arg); read(arg, *) gen
  call get_command_argument(2, arg); read(arg, *) seed
  call get_command_argument(3, arg); read(arg, *) single
  call get_command_argument(4, arg); read(arg, *) l2max
  call get_command_argument(5, outfile)
  rng_single = (single /= 0)
  nmax = 2_8**l2max

  if (single /= 0) then
     nbits = 1065353216_8          ! bit patterns of floats in [0,1)
  else
     select case (gen)
     case (G_RANDU, G_RANDM, G_MINSTD); nbits = 2147483648_8
     case (G_RAND);                     nbits = 4194304_8
     case default; stop 'generator not supported by the bitmap test'
     end select
  end if
  allocate(bm(0:(nbits - 1) / 64))
  bm = 0_8
  call rng_init(gen, seed)
  open(20, file=outfile, status='replace')
  write(20, '(a,a,a,i0,a,i0)') '# ', trim(gen_name(gen)), ' seed ', seed, ' single ', single
  write(20, '(a)') '# N D'
  d = 0
  icp = 40                     ! checkpoints at N = 2^(icp/4)
  nextcp = nint(2.0d0**(icp / 4.0d0), 8)
  do n = 1, nmax
     u = rng_next(gen)
     if (single /= 0) then
        us = real(u, 4)
        key = int(transfer(us, 1_4), 8)
     else
        select case (gen)
        case (G_RAND); key = nint(u * 4194304.0d0, 8)
        case default;  key = rng_state_int(gen)
        end select
     end if
     w = ishft(key, -6)
     ib = int(iand(key, 63_8))
     if (btest(bm(w), ib)) then
        d = d + 1
     else
        bm(w) = ibset(bm(w), ib)
     end if
     if (n == nextcp) then
        write(20, '(i0,1x,i0)') n, d
        do while (nextcp <= n)
           icp = icp + 1
           nextcp = nint(2.0d0**(icp / 4.0d0), 8)
        end do
     end if
  end do
  close(20)
end program spacing_dups
