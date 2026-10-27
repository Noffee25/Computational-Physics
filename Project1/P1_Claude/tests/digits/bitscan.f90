! =====================================================================
!  bitscan.f90 : bit-level statistics of the first 53 binary digits
!
!  usage: bitscan gid prec N seed outfile dumpfile ndump
!
!  b_j(n) = floor(x_n 2^j) mod 2, j = 1..53, read exactly from the
!  integer m = x 2^53 (exact because multiplication by 2^53 is exact).
!  Bit j of x is integer bit 53-j of m.
!  Per bit j it accumulates
!    ones      : number of n with b_j(n) = 1
!    a11(L)    : number of n with b_j(n) = b_j(n-L) = 1, L = 1, 2, 4
!    v250      : number of n with b_j(n) XOR b_j(n-103) XOR b_j(n-250) = 1
!                (R250 recurrence identity)
!    vtent     : number of n with b_j(n) XOR b_{j+1}(n-1) XOR b_1(n-1) = 1
!                (tent-map / logistic shift identity), j = 1..52
!    vcopy     : number of n with b_{j+31}(n) XOR b_j(n) = 1, j = 1..22
!                (period-31 binary expansion of k/(2^31-1))
!  and dumps the first ndump integers m to dumpfile (bit-plane images).
! =====================================================================
program bitscan
  use rnglib
  implicit none
  integer(8), parameter :: TWO53 = 9007199254740992_8, MASK53 = TWO53 - 1_8
  integer(8), parameter :: MASK22 = 4194303_8
  integer :: gid, iprec, b, ndump, il
  integer(8) :: nd, seed, n, m, mp, v, ring(0:255)
  integer(8) :: ones(0:52), a11(0:52, 3), v250(0:52), vtent(0:52), vcopy(0:52)
  integer(8) :: n250, nlag(3)
  integer, parameter :: lags(3) = [1, 2, 4]
  real(8) :: x
  character(len=256) :: arg, outf, dumpf

  call get_command_argument(1, arg); read(arg, *) gid
  call get_command_argument(2, arg); read(arg, *) iprec
  call get_command_argument(3, arg); read(arg, *) nd
  call get_command_argument(4, arg); read(arg, *) seed
  call get_command_argument(5, outf)
  call get_command_argument(6, dumpf)
  call get_command_argument(7, arg); read(arg, *) ndump

  rng_single = (iprec == 1)
  call rng_init(gid, seed)
  ones = 0; a11 = 0; v250 = 0; vtent = 0; vcopy = 0; n250 = 0; nlag = 0
  ring = 0
  open(11, file=trim(dumpf), status='replace')
  do n = 0, nd - 1
     x = rng_next(gid)
     m = int(x * 9007199254740992.0d0, 8)
     if (n < ndump) write(11, '(i0)') m
     do b = 0, 52
        ones(b) = ones(b) + ibits(m, b, 1)
     end do
     do il = 1, 3
        if (n >= lags(il)) then
           v = iand(m, ring(iand(n - lags(il), 255_8)))
           nlag(il) = nlag(il) + 1
           do b = 0, 52
              a11(b, il) = a11(b, il) + ibits(v, b, 1)
           end do
        end if
     end do
     if (n >= 1) then
        mp = ring(iand(n - 1, 255_8))
        v = ieor(m, iand(shiftl(mp, 1), MASK53))
        if (btest(mp, 52)) v = ieor(v, MASK53)
        do b = 1, 52
           vtent(b) = vtent(b) + ibits(v, b, 1)
        end do
     end if
     if (n >= 250) then
        v = ieor(ieor(m, ring(iand(n - 103, 255_8))), ring(iand(n - 250, 255_8)))
        n250 = n250 + 1
        do b = 0, 52
           v250(b) = v250(b) + ibits(v, b, 1)
        end do
     end if
     v = iand(ieor(m, shiftr(m, 31)), MASK22)
     do b = 0, 21
        vcopy(b) = vcopy(b) + ibits(v, b, 1)
     end do
     ring(iand(n, 255_8)) = m
  end do
  close(11)

  open(10, file=trim(outf), status='replace')
  write(10, '(a,i3,a,i2,a,i12,a,i12,a,3i12,a,i12)') '# gid ', gid, ' prec ', iprec, ' N ', nd, &
       ' seed ', seed, ' nlag ', nlag, ' n250 ', n250
  write(10, '(a)') '# j ones a11_L1 a11_L2 a11_L4 v250 vtent vcopy'
  do b = 52, 0, -1          ! bit j = 53 - b
     write(10, '(i3,7i12)') 53 - b, ones(b), a11(b, 1), a11(b, 2), a11(b, 3), v250(b), &
          vtent(b), vcopy(b)
  end do
  close(10)
  write(*, '(a,a,a,i2)') 'done ', trim(gen_name(gid)), ' prec ', iprec
end program bitscan
