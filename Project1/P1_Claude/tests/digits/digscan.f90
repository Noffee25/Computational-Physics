! =====================================================================
!  digscan.f90 : decimal digit frequencies at positions p = 1..PMAX
!
!  usage: digscan gid prec N seed outfile
!     gid  : generator id (rnglib)
!     prec : 0 = double precision output, 1 = rng_single (REAL*4 values)
!
!  For every draw x the program computes, at every position p,
!   STO : the EXACT p-th decimal digit of the stored binary number x
!         (x = m 2^-B is a dyadic rational; digits by exact integer
!          arithmetic r <- 10 r, d = r >> B, r <- r mod 2^B)
!   INT : the exact digit of the value the generator intends
!         (k/(2^31-1) for MINSTD, k/10^8 for MIDSQ, from the integer
!          state; identical to STO for every other generator)
!   NAI : the naive digit  mod(int(x*10^p), 10)  in the working
!         precision (double: p <= 18; single: REAL*4 product, p <= 12)
!  and accumulates the 10 x PMAX count tables plus mismatch counts
!  (NAI vs STO, STO vs INT).
! =====================================================================
program digscan
  use rnglib
  implicit none
  integer, parameter :: PMAX = 60, PND = 18, PNS = 12
  integer :: gid, iprec, p, dn, i
  integer(8) :: n, nd, seed, k, r, q, modl
  integer(8) :: cnt_int(0:9, PMAX), cnt_sto(0:9, PMAX), cnt_nai(0:9, PND)
  integer(8) :: mis_nai(PND), mis_int(PMAX)
  integer :: dig_int(PMAX), dig_sto(PMAX)
  real(8) :: x, p10(0:PND), t0, t1
  real(4) :: x4, y4, p10s(0:PNS)
  character(len=256) :: arg, outf
  logical :: has_int

  call get_command_argument(1, arg); read(arg, *) gid
  call get_command_argument(2, arg); read(arg, *) iprec
  call get_command_argument(3, arg); read(arg, *) nd
  call get_command_argument(4, arg); read(arg, *) seed
  call get_command_argument(5, outf)

  rng_single = (iprec == 1)
  call rng_init(gid, seed)
  has_int = (iprec == 0) .and. (gid == G_MINSTD .or. gid == G_MIDSQ)
  if (gid == G_MINSTD) modl = 2147483647_8
  if (gid == G_MIDSQ)  modl = 100000000_8

  p10(0) = 1.0d0
  do p = 1, PND
     p10(p) = p10(p-1) * 10.0d0          ! exact: 10^p < 2^53 for p <= 22
  end do
  p10s(0) = 1.0
  do p = 1, PNS
     p10s(p) = real(p10(p), 4)           ! REAL*4 10^p (exact up to 10^10)
  end do

  cnt_int = 0; cnt_sto = 0; cnt_nai = 0; mis_nai = 0; mis_int = 0
  call cpu_time(t0)
  do n = 1, nd
     x = rng_next(gid)
     call exact_digits(x, dig_sto, PMAX)
     do p = 1, PMAX
        cnt_sto(dig_sto(p), p) = cnt_sto(dig_sto(p), p) + 1
     end do
     if (has_int) then
        r = rng_state_int(gid)
        do p = 1, PMAX
           r = r * 10_8
           q = r / modl
           r = r - q * modl
           dig_int(p) = int(q)
           cnt_int(dig_int(p), p) = cnt_int(dig_int(p), p) + 1
           if (dig_int(p) /= dig_sto(p)) mis_int(p) = mis_int(p) + 1
        end do
     end if
     if (iprec == 0) then
        do p = 1, PND
           dn = int(mod(int(x * p10(p), 8), 10_8))
           cnt_nai(dn, p) = cnt_nai(dn, p) + 1
           if (dn /= dig_sto(p)) mis_nai(p) = mis_nai(p) + 1
        end do
     else
        x4 = real(x, 4)
        do p = 1, PNS
           y4 = x4 * p10s(p)                 ! single precision product
           dn = int(mod(int(y4, 8), 10_8))
           cnt_nai(dn, p) = cnt_nai(dn, p) + 1
           if (dn /= dig_sto(p)) mis_nai(p) = mis_nai(p) + 1
        end do
     end if
  end do
  call cpu_time(t1)
  if (.not. has_int) cnt_int = cnt_sto

  open(10, file=trim(outf), status='replace')
  write(10, '(a,i3,a,i2,a,i12,a,i12,a,f8.2)') '# gid ', gid, ' prec ', iprec, ' N ', nd, &
       ' seed ', seed, ' cpu ', t1 - t0
  do p = 1, PMAX
     write(10, '(a,i3,10i12)') 'STO ', p, (cnt_sto(i, p), i = 0, 9)
  end do
  do p = 1, PMAX
     write(10, '(a,i3,10i12)') 'INT ', p, (cnt_int(i, p), i = 0, 9)
  end do
  if (iprec == 0) then
     do p = 1, PND
        write(10, '(a,i3,10i12)') 'NAI ', p, (cnt_nai(i, p), i = 0, 9)
     end do
     do p = 1, PND
        write(10, '(a,i3,i12)') 'MISNAI ', p, mis_nai(p)
     end do
  else
     do p = 1, PNS
        write(10, '(a,i3,10i12)') 'NAI ', p, (cnt_nai(i, p), i = 0, 9)
     end do
     do p = 1, PNS
        write(10, '(a,i3,i12)') 'MISNAI ', p, mis_nai(p)
     end do
  end if
  do p = 1, PMAX
     write(10, '(a,i3,i12)') 'MISINT ', p, mis_int(p)
  end do
  close(10)
  write(*, '(a,a,a,i2,a,f8.2,a)') 'done ', trim(gen_name(gid)), ' prec ', iprec, ' cpu ', t1 - t0, ' s'

contains

  ! exact decimal digits 1..np of a binary floating point number 0 <= x < 1
  subroutine exact_digits(x, dig, np)
    real(8), intent(in) :: x
    integer, intent(in) :: np
    integer, intent(out) :: dig(np)
    integer :: e, b, t, p
    integer(8) :: m, r, mask
    integer(16) :: r16, mask16
    if (x <= 0.0d0) then
       dig = 0
       return
    end if
    e = exponent(x)                         ! x = f 2^e, 0.5 <= f < 1
    m = int(scale(fraction(x), 53), 8)      ! exact 53-bit integer
    b = 53 - e                              ! x = m / 2^b exactly
    t = trailz(m)
    m = shiftr(m, t)
    b = b - t
    if (b <= 59) then
       mask = shiftl(1_8, b) - 1_8
       r = m
       do p = 1, np
          r = r * 10_8
          dig(p) = int(shiftr(r, b))
          r = iand(r, mask)
       end do
    else
       if (b > 122) then                    ! x < 2^-69: never happens in practice
          m = shiftr(m, b - 122)
          b = 122
       end if
       mask16 = shiftl(1_16, b) - 1_16
       r16 = int(m, 16)
       do p = 1, np
          r16 = r16 * 10_16
          dig(p) = int(shiftr(r16, b))
          r16 = iand(r16, mask16)
       end do
    end if
  end subroutine exact_digits

end program digscan
