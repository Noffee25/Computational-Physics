! =====================================================================
!  bitperiod.f90 : exact periods of the individual bits and decimal
!                  digits of the power-of-two-modulus LCGs
!
!  usage: bitperiod gid seed outfile        (gid = 1 RANDU, 2 RANDM, 4 MINSTD)
!
!  A. First return times R(i) of the state residue x_n mod 2^i to
!     x_0 mod 2^i, i = 1..31, over the full orbit (run until x_n = x_0),
!     and a check that at n = R(i)/2 exactly bit i-1 differs.
!  B. Brute-force minimal period of each bit sequence b_i(n) (bit i of
!     the integer state, LSB = 0) for i = 0..22, and of each exact
!     decimal digit sequence d_p(n) = floor(10^p k_n / 2^31) mod 10,
!     p = 8..30, using L = 2^23 consecutive states.
!  C. Half-period identity of the digits: fraction of n with
!     d_p(n + P/2) = d_p(n) + 5 (mod 10).
!  D. Lagged autocorrelation of bits at lags 1, 2, 4, ..., 2^12.
! =====================================================================
program bitperiod
  use rnglib
  implicit none
  integer, parameter :: LL = 8388608           ! 2^23 stored states
  integer(8) :: seed, x0, x, n, v, rr(31), period, smax
  integer :: gid, i, nxt, t, p, il, b, ok
  integer(4), allocatable :: st(:)
  integer(1), allocatable :: sq(:)
  integer(8) :: per_bit(0:22), per_dig(8:30), nmatch, nflip(31)
  integer(16) :: k16, p10
  real(8) :: x8, f1, f0, cor, ones, a11
  character(len=256) :: arg, outf

  call get_command_argument(1, arg); read(arg, *) gid
  call get_command_argument(2, arg); read(arg, *) seed
  call get_command_argument(3, outf)
  open(10, file=trim(outf), status='replace')
  write(10, '(a,a,a,i12)') '# generator ', trim(gen_name(gid)), ' seed ', seed

  ! ---------------- A. residue return times over the whole orbit
  if (gid /= G_MINSTD) then
     call rng_init(gid, seed)
     x0 = modulo(seed, 2147483648_8)
     rr = 0; nxt = 1; n = 0; period = -1
     do
        x8 = rng_next(gid); x = rng_state_int(gid); n = n + 1
        v = ieor(x, x0)
        if (v == 0) then
           period = n
           do while (nxt <= 31)
              rr(nxt) = n; nxt = nxt + 1
           end do
           exit
        end if
        t = trailz(v)
        do while (nxt <= t)
           rr(nxt) = n; nxt = nxt + 1
        end do
        if (n > 2147483648_8) exit
     end do
     ! check flip of bit i-1 at half of R(i)
     call rng_init(gid, seed)
     nflip = 0
     do n = 1, period / 2
        x8 = rng_next(gid); x = rng_state_int(gid)
        do i = 2, 31
           if (n == rr(i) / 2) then
              v = ieor(x, x0)
              if (trailz(v) == i - 1) nflip(i) = 1
           end if
        end do
     end do
     write(10, '(a,i14)') 'FULLPERIOD ', period
     do i = 1, 31
        write(10, '(a,i3,i14,i3)') 'RETURN ', i, rr(i), nflip(i)
     end do
  end if

  ! ---------------- B. brute-force minimal periods from 2^23 states
  allocate(st(0:LL-1), sq(0:LL-1))
  call rng_init(gid, seed)
  do n = 0, LL - 1
     x8 = rng_next(gid)
     st(n) = int(rng_state_int(gid), 4)
  end do
  smax = LL / 4
  do b = 0, 22
     do n = 0, LL - 1
        sq(n) = int(ibits(st(n), b, 1), 1)
     end do
     per_bit(b) = minper(sq, LL, smax)
     write(10, '(a,i3,i14)') 'BITPERIOD ', b, per_bit(b)
  end do
  if (gid /= G_MINSTD) then
     do p = 8, 30
        p10 = 5_16**p                       ! 10^p / 2^31 = 5^p / 2^(31-p)
        do n = 0, LL - 1
           k16 = int(st(n), 16) * p10
           sq(n) = int(mod(shiftr(k16, 31 - p), 10_16), 1)
        end do
        per_dig(p) = minper(sq, LL, smax)
        ! C. half-period identity
        nmatch = 0
        if (per_dig(p) > 1) then
           do n = 0, LL - 1 - per_dig(p) / 2
              if (mod(int(sq(n + per_dig(p) / 2)) - int(sq(n)) + 10, 10) == 5) nmatch = nmatch + 1
           end do
        end if
        write(10, '(a,i3,i14,i14,i14)') 'DIGPERIOD ', p, per_dig(p), nmatch, LL - per_dig(p) / 2
     end do
  end if

  ! ---------------- D. lagged autocorrelation of each state bit
  do b = 0, 30
     do n = 0, LL - 1
        sq(n) = int(ibits(st(n), b, 1), 1)
     end do
     ones = dble(count(sq == 1)) / dble(LL)
     do il = 0, 12
        a11 = 0.0d0
        do n = 0, LL - 1 - 2**il
           a11 = a11 + dble(iand(sq(n), sq(n + 2**il)))
        end do
        a11 = a11 / dble(LL - 2**il)
        if (ones > 0.0d0 .and. ones < 1.0d0) then
           cor = (a11 - ones * ones) / (ones * (1.0d0 - ones))
        else
           cor = 2.0d0                      ! flag: constant bit
        end if
        write(10, '(a,i3,i3,f12.6)') 'LAGCOR ', b, il, cor
     end do
  end do
  close(10)
  write(*, '(a,a)') 'done ', trim(gen_name(gid))

contains

  ! smallest s in 1..smax with q(n+s) = q(n) for all valid n; -1 if none
  integer(8) function minper(q, nq, smax)
    integer, intent(in) :: nq
    integer(1), intent(in) :: q(0:nq-1)
    integer(8), intent(in) :: smax
    integer(8) :: s, j
    logical :: good
    minper = -1
    do s = 1, smax
       good = .true.
       do j = 0, nq - 1 - s
          if (q(j + s) /= q(j)) then
             good = .false.
             exit
          end if
       end do
       if (good) then
          minper = s
          return
       end if
    end do
  end function minper

end program bitperiod
