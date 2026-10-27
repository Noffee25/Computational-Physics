! =====================================================================
!  bm.f90 : linear complexity profile of a single bit stream
!           (Berlekamp-Massey algorithm over GF(2))
!
!  usage: bm gid seed jbit L outfile
!
!  The stream is s_n = b_j(x_n) = floor(x_n 2^j) mod 2, n = 0..L-1.
!  The linear complexity LC(n) is the length of the shortest linear
!  feedback shift register that generates s_0..s_{n-1}.  For a random
!  sequence LC(n) follows n/2 closely (deviations of a few bits); for a
!  sequence produced by a GF(2)-linear recurrence of order k it stops
!  at k once n > 2k.  Writes LC(n) every 10 steps and the final value.
! =====================================================================
program bm
  use rnglib
  implicit none
  integer :: gid, jbit, lcur, mlast, ldb, i, nn
  integer(8) :: seed, m
  integer :: ll
  integer(1), allocatable :: s(:), sr(:), c(:), bb(:), tt(:)
  integer(1) :: d
  real(8) :: x
  character(len=256) :: arg, outf
  call get_command_argument(1, arg); read(arg, *) gid
  call get_command_argument(2, arg); read(arg, *) seed
  call get_command_argument(3, arg); read(arg, *) jbit
  call get_command_argument(4, arg); read(arg, *) ll
  call get_command_argument(5, outf)
  call rng_init(gid, seed)
  allocate(s(0:ll-1), sr(0:ll-1), c(0:ll), bb(0:ll), tt(0:ll))
  do nn = 0, ll - 1
     x = rng_next(gid)
     m = int(x * 9007199254740992.0d0, 8)
     s(nn) = int(ibits(m, 53 - jbit, 1), 1)
  end do
  do nn = 0, ll - 1
     sr(nn) = s(ll - 1 - nn)             ! reversed copy: s(n-i) = sr(ll-1-n+i)
  end do
  c = 0; bb = 0; c(0) = 1; bb(0) = 1
  lcur = 0; mlast = -1; ldb = 0
  open(10, file=trim(outf), status='replace')
  write(10, '(a,a,a,i12,a,i3,a,i8)') '# ', trim(gen_name(gid)), ' seed ', seed, ' bit ', jbit, ' L ', ll
  do nn = 0, ll - 1
     d = s(nn)
     do i = 1, lcur
        d = ieor(d, iand(c(i), sr(ll - 1 - nn + i)))
     end do
     if (d /= 0) then
        if (2 * lcur <= nn) then
           tt(0:lcur) = c(0:lcur)
           do i = 0, ldb
              c(nn - mlast + i) = ieor(c(nn - mlast + i), bb(i))
           end do
           bb(0:lcur) = tt(0:lcur)
           ldb = lcur
           lcur = nn + 1 - lcur
           mlast = nn
        else
           do i = 0, ldb
              c(nn - mlast + i) = ieor(c(nn - mlast + i), bb(i))
           end do
        end if
     end if
     if (mod(nn + 1, 10) == 0) write(10, '(i8,i8)') nn + 1, lcur
  end do
  write(10, '(a,i8)') 'FINAL ', lcur
  close(10)
  write(*, '(a,a,a,i3,a,i8)') 'done ', trim(gen_name(gid)), ' bit ', jbit, ' LC ', lcur
end program bm
