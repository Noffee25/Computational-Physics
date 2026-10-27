! =====================================================================
!  rnglib.f90 : shared random number generator library for Project 1
!
!  Every test program in P1_Claude uses this module so that all results
!  come from one, verified implementation of each generator.
!
!  Usage (free or fixed form):
!      use rnglib
!      call rng_init(G_RANDU, 777_8)
!      x = rng_next(G_RANDU)          ! real(8), 0 <= x < 1
!
!  Set rng_single = .true. to round every output to single precision
!  (REAL*4) before it is returned, emulating a REAL variable holding
!  the random number.
!
!  Compile:  gfortran -O2 -c rnglib.f90   (then link rnglib.o)
! =====================================================================
module rnglib
  implicit none
  private

  integer, parameter, public :: G_RANDU  = 1   ! IBM RANDU, a=65539, m=2^31
  integer, parameter, public :: G_RANDM  = 2   ! course RANDM, a=6539, m=2^31
  integer, parameter, public :: G_RAND   = 3   ! gfortran intrinsic RAND()
  integer, parameter, public :: G_MINSTD = 4   ! Park-Miller, a=16807, m=2^31-1
  integer, parameter, public :: G_R250   = 5   ! Kirkpatrick-Stoll R250 (XOR, lags 103/250)
  integer, parameter, public :: G_RANNUM = 6   ! gfortran RANDOM_NUMBER (xoshiro256**)
  integer, parameter, public :: G_MT     = 7   ! Mersenne Twister MT19937 (53-bit output)
  integer, parameter, public :: G_MIDSQ  = 8   ! von Neumann middle-square, 8 digits
  integer, parameter, public :: G_LOGIS  = 9   ! logistic map x->4x(1-x), mapped to uniform
  integer, parameter, public :: NGEN     = 9

  character(len=8), parameter, public :: gen_name(NGEN) = [character(len=8) :: &
       'RANDU', 'RANDM', 'RAND', 'MINSTD', 'R250', 'RANNUM', 'MT19937', &
       'MIDSQ', 'LOGISTIC']

  logical, public :: rng_single = .false.

  public :: rng_init, rng_next, rng_fill, rng_state_int, mt_u32

  integer(8), parameter :: M31  = 2147483648_8      ! 2^31
  integer(8), parameter :: PM_M = 2147483647_8      ! 2^31 - 1
  integer(8), parameter :: MASK32 = 4294967295_8    ! 2^32 - 1
  real(8),    parameter :: TWO32 = 4294967296.0d0
  real(8),    parameter :: ONE_MINUS = 1.0d0 - epsilon(1.0d0)
  real(8),    parameter :: ONE_MINUS_SP = 0.99999994039535522d0  ! largest REAL*4 below 1

  ! --- LCG states
  integer(8) :: st_randu = 777_8, st_randm = 777_8, st_minstd = 1_8, st_midsq = 12345678_8
  ! --- logistic map state
  real(8) :: st_logis = 0.3141592653589793d0
  ! --- R250 state
  integer(8) :: r250_buf(0:249)
  integer    :: r250_idx = 0
  logical    :: r250_ready = .false.
  ! --- MT19937 state
  integer(8) :: mt(0:623)
  integer    :: mti = 625

contains

  ! -------------------------------------------------------------------
  subroutine rng_init(id, seed)
    integer,    intent(in) :: id
    integer(8), intent(in) :: seed
    select case (id)
    case (G_RANDU);  st_randu  = modulo(seed, M31)
    case (G_RANDM);  st_randm  = modulo(seed, M31)
    case (G_RAND);   call srand(int(seed, 4))
    case (G_MINSTD)
       st_minstd = modulo(seed, PM_M)
       if (st_minstd == 0) st_minstd = 1
    case (G_R250);   call r250_init(seed)
    case (G_RANNUM); call rannum_init(seed)
    case (G_MT);     call mt_init(iand(seed, MASK32))
    case (G_MIDSQ);  st_midsq = modulo(seed, 100000000_8)
    case (G_LOGIS)
       st_logis = modulo(dble(seed) * 0.6180339887498949d0, 1.0d0)
       if (st_logis <= 0.0d0 .or. st_logis >= 1.0d0) st_logis = 0.3141592653589793d0
    case default
       stop 'rng_init: unknown generator id'
    end select
  end subroutine rng_init

  ! -------------------------------------------------------------------
  function rng_next(id) result(x)
    integer, intent(in) :: id
    real(8) :: x
    real(4) :: xs
    select case (id)
    case (G_RANDU)
       st_randu = modulo(65539_8 * st_randu, M31)
       x = dble(st_randu) / dble(M31)
    case (G_RANDM)
       st_randm = modulo(6539_8 * st_randm, M31)
       x = dble(st_randm) / dble(M31)
    case (G_RAND)
       xs = rand()
       x = dble(xs)
    case (G_MINSTD)
       st_minstd = modulo(16807_8 * st_minstd, PM_M)
       x = dble(st_minstd) / dble(PM_M)
    case (G_R250)
       x = dble(r250_word()) / TWO32
    case (G_RANNUM)
       call random_number(x)
    case (G_MT)
       x = mt_res53()
    case (G_MIDSQ)
       st_midsq = modulo((st_midsq * st_midsq) / 10000_8, 100000000_8)
       x = dble(st_midsq) / 1.0d8
    case (G_LOGIS)
       st_logis = 4.0d0 * st_logis * (1.0d0 - st_logis)
       x = (2.0d0 / 3.141592653589793d0) * asin(sqrt(st_logis))
    case default
       stop 'rng_next: unknown generator id'
    end select
    if (x >= 1.0d0) x = ONE_MINUS
    if (x < 0.0d0)  x = 0.0d0
    if (rng_single) then
       x = dble(real(x, 4))
       if (x >= 1.0d0) x = ONE_MINUS_SP
    end if
  end function rng_next

  ! -------------------------------------------------------------------
  subroutine rng_fill(id, x, n)
    integer, intent(in)  :: id, n
    real(8), intent(out) :: x(n)
    integer :: i
    do i = 1, n
       x(i) = rng_next(id)
    end do
  end subroutine rng_fill

  ! -------------------------------------------------------------------
  ! Internal integer state of the LCG-type generators (RANDU, RANDM,
  ! MINSTD, MIDSQ). Returns -1 for generators without a single integer
  ! state. Useful for exact period measurements.
  function rng_state_int(id) result(s)
    integer, intent(in) :: id
    integer(8) :: s
    select case (id)
    case (G_RANDU);  s = st_randu
    case (G_RANDM);  s = st_randm
    case (G_MINSTD); s = st_minstd
    case (G_MIDSQ);  s = st_midsq
    case default;    s = -1_8
    end select
  end function rng_state_int

  ! =================== R250 ==========================================
  ! x_n = x_{n-103} XOR x_{n-250} on 32-bit words (Kirkpatrick & Stoll
  ! 1981). Seeded as in W. L. Maier's r250.c: buffer filled from a
  ! Park-Miller LCG, then a diagonal of bits forced on to guarantee the
  ! 32 bit columns are linearly independent.
  subroutine r250_init(seed)
    integer(8), intent(in) :: seed
    integer(8) :: s, msb, mask
    integer :: j, k
    s = modulo(seed, PM_M); if (s == 0) s = 1
    do j = 0, 249
       s = modulo(16807_8 * s, PM_M)
       r250_buf(j) = s                       ! 31 random bits
    end do
    do j = 0, 249
       s = modulo(16807_8 * s, PM_M)
       if (s > 1073741823_8) r250_buf(j) = ior(r250_buf(j), 2147483648_8)  ! set bit 31
    end do
    msb  = 2147483648_8
    mask = MASK32
    do j = 0, 31
       k = 7 * j + 3
       r250_buf(k) = ior(iand(r250_buf(k), mask), msb)
       mask = ishft(mask, -1)
       msb  = ishft(msb, -1)
    end do
    r250_idx = 0
    r250_ready = .true.
  end subroutine r250_init

  function r250_word() result(w)
    integer(8) :: w
    integer :: j
    if (.not. r250_ready) call r250_init(1_8)
    ! position r250_idx holds x_{n-250}; position r250_idx-103 holds x_{n-103}
    j = r250_idx - 103
    if (j < 0) j = j + 250
    w = ieor(r250_buf(r250_idx), r250_buf(j))
    r250_buf(r250_idx) = w
    r250_idx = r250_idx + 1
    if (r250_idx >= 250) r250_idx = 0
  end function r250_word

  ! =================== RANDOM_NUMBER seeding ==========================
  subroutine rannum_init(seed)
    integer(8), intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    integer(8) :: z
    call random_seed(size=n)
    allocate(put(n))
    z = seed
    do i = 1, n
       z = splitmix_step(z)
       put(i) = int(iand(z, MASK32) - 2147483648_8, 4)
       if (put(i) == 0) put(i) = 1
    end do
    call random_seed(put=put)
    deallocate(put)
  end subroutine rannum_init

  ! A 31-bit mixing step (pure integer arithmetic, no overflow) used only
  ! to spread a user seed over RANDOM_NUMBER's seed array.
  function splitmix_step(z) result(r)
    integer(8), intent(in) :: z
    integer(8) :: r
    r = modulo(z * 1103515245_8 + 12345_8, M31)
    r = modulo(ieor(r, ishft(r, -13)) * 69069_8 + 1_8, M31)
    r = ior(ishft(r, 1), iand(ishft(r, -7), 1_8))
  end function splitmix_step

  ! =================== MT19937 =======================================
  subroutine mt_init(s)
    integer(8), intent(in) :: s
    integer :: i
    mt(0) = iand(s, MASK32)
    do i = 1, 623
       mt(i) = iand(1812433253_8 * ieor(mt(i-1), ishft(mt(i-1), -30)) + int(i, 8), MASK32)
    end do
    mti = 624
  end subroutine mt_init

  function mt_u32() result(y)
    integer(8) :: y
    integer :: kk
    integer(8), parameter :: UPPER = 2147483648_8, LOWER = 2147483647_8
    integer(8), parameter :: MATA = 2567483615_8          ! 0x9908b0df
    integer(8), parameter :: TB = 2636928640_8, TC = 4022730752_8  ! 0x9d2c5680, 0xefc60000
    if (mti >= 624) then
       if (mti == 625) call mt_init(5489_8)
       do kk = 0, 623
          y = ior(iand(mt(kk), UPPER), iand(mt(mod(kk+1, 624)), LOWER))
          mt(kk) = ieor(mt(mod(kk+397, 624)), ishft(y, -1))
          if (btest(y, 0)) mt(kk) = ieor(mt(kk), MATA)
       end do
       mti = 0
    end if
    y = mt(mti)
    mti = mti + 1
    y = ieor(y, ishft(y, -11))
    y = ieor(y, iand(ishft(y, 7),  TB))
    y = ieor(y, iand(ishft(y, 15), TC))
    y = ieor(y, ishft(y, -18))
    y = iand(y, MASK32)
  end function mt_u32

  function mt_res53() result(x)
    real(8) :: x
    integer(8) :: a, b
    a = ishft(mt_u32(), -5)
    b = ishft(mt_u32(), -6)
    x = (dble(a) * 67108864.0d0 + dble(b)) / 9007199254740992.0d0
  end function mt_res53

end module rnglib
