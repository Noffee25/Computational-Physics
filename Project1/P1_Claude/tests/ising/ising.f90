! =====================================================================
!  ising.f90 : 2D Ising model on an L x L torus, Wolff single-cluster
!  and Metropolis updates, driven by the generators of rnglib plus a
!  family of generalized feedback shift register (GFSR) variants.
!
!  Reads a namelist &run from stdin. Writes one output file with the
!  per-batch integer sums (exact, no roundoff) and the energy
!  autocorrelation sums. All analysis is done in Python (analyze.py).
!
!  Sources (src):
!    1..9  rnglib generators (RANDU ... LOGISTIC)
!    20    own GFSR/two- or four-tap XOR shift register, lags(1:ntap)
!    21    same, decimated: keep every dec-th word
!    22    XOR of two independent GFSR streams (seeds seed, seed2)
!    23    sum mod 2^32 of two independent GFSR streams
!    24    additive lagged Fibonacci x_n = x_{n-p} + x_{n-q} mod 2^32
!    25    GFSR XOR MT19937 32-bit word
! =====================================================================
module rsrc
  use rnglib
  implicit none
  integer, parameter :: S_GFSR = 20, S_DEC = 21, S_XOR2 = 22, S_ADD2 = 23, S_ALFG = 24, S_GMT = 25
  real(8), parameter :: TWO32 = 4294967296.0d0
  integer(8), parameter :: MASK32 = 4294967295_8
  integer :: src = 7
  integer :: dec = 1
  integer(8) :: ndraw = 0

  type gfsr_t
     integer :: q = 0, ntap = 0, idx = 0
     integer :: off(3) = 0
     logical :: additive = .false.
     integer(8), allocatable :: buf(:)
  end type gfsr_t
  type(gfsr_t) :: g1, g2

contains

  ! Seeding identical to rnglib's r250_init (W. L. Maier's r250.c scheme)
  subroutine gfsr_init(g, lags, ntap, seed, additive)
    type(gfsr_t), intent(inout) :: g
    integer, intent(in) :: lags(4), ntap
    integer(8), intent(in) :: seed
    logical, intent(in) :: additive
    integer(8) :: s, msb, mask
    integer :: j, k
    g%q = lags(ntap); g%ntap = ntap; g%additive = additive
    if (allocated(g%buf)) deallocate(g%buf)
    allocate(g%buf(0:g%q-1))
    g%off = 0
    do k = 1, ntap - 1
       g%off(k) = g%q - lags(k)
    end do
    s = modulo(seed, 2147483647_8); if (s == 0) s = 1
    do j = 0, g%q - 1
       s = modulo(16807_8 * s, 2147483647_8); g%buf(j) = s
    end do
    do j = 0, g%q - 1
       s = modulo(16807_8 * s, 2147483647_8)
       if (s > 1073741823_8) g%buf(j) = ior(g%buf(j), 2147483648_8)
    end do
    msb = 2147483648_8; mask = MASK32
    do j = 0, 31
       k = 7 * j + 3
       if (k >= g%q) exit           ! short registers (q <= 220): no diagonal
       g%buf(k) = ior(iand(g%buf(k), mask), msb)
       mask = ishft(mask, -1); msb = ishft(msb, -1)
    end do
    g%idx = 0
  end subroutine gfsr_init

  function gfsr_word(g) result(w)
    type(gfsr_t), intent(inout) :: g
    integer(8) :: w
    integer :: k, j
    w = g%buf(g%idx)
    do k = 1, g%ntap - 1
       j = g%idx + g%off(k); if (j >= g%q) j = j - g%q
       if (g%additive) then
          w = iand(w + g%buf(j), MASK32)
       else
          w = ieor(w, g%buf(j))
       end if
    end do
    g%buf(g%idx) = w
    g%idx = g%idx + 1; if (g%idx >= g%q) g%idx = 0
  end function gfsr_word

  function u() result(x)
    real(8) :: x
    integer(8) :: w
    integer :: k
    ndraw = ndraw + 1
    select case (src)
    case (1:9)
       x = rng_next(src)
    case (S_GFSR, S_ALFG)
       x = dble(gfsr_word(g1)) / TWO32
    case (S_DEC)
       do k = 1, dec
          w = gfsr_word(g1)
       end do
       x = dble(w) / TWO32
    case (S_XOR2)
       x = dble(ieor(gfsr_word(g1), gfsr_word(g2))) / TWO32
    case (S_ADD2)
       x = dble(iand(gfsr_word(g1) + gfsr_word(g2), MASK32)) / TWO32
    case (S_GMT)
       x = dble(ieor(gfsr_word(g1), mt_u32())) / TWO32
    case default
       stop 'u: bad src'
    end select
  end function u

end module rsrc

! =====================================================================
program ising
  use rnglib
  use rsrc
  implicit none
  character(len=8)   :: alg = 'wolff'      ! 'wolff' | 'metro' | 'mrand' | 'mskip' | 'check' | 'diag'
  character(len=4)   :: order = 'lifo'     ! Wolff cluster growth: 'lifo' (stack) | 'fifo' (queue)
  character(len=200) :: outfile = 'out.dat'
  integer    :: L = 16, ntap = 2, nbatch = 1000, lags(4) = [103, 250, 0, 0]
  integer(8) :: ntherm = 10000, nmeas = 1000000, seed = 777, seed2 = 1234567
  real(8)    :: K = 0.44068679350977151d0
  logical    :: ordered_start = .false.
  namelist /run/ alg, order, outfile, L, ntap, nbatch, lags, ntherm, nmeas, seed, seed2, K, src, dec, ordered_start

  integer, parameter :: NACF = 64
  integer :: N, i, ib, nper, csize
  integer, allocatable :: sp(:), nbr(:,:), stk(:)
  integer(8) :: t, B, M, sB, sB2, sC, sAM, sM2, acf(0:NACF), ring(0:NACF-1), it
  real(8) :: padd, w4, w8, t0, t1
  integer :: ix, iy
  ! diag-mode state
  integer, parameter :: LMAX = 400
  integer :: tsite(0:LMAX), tclus(0:LMAX), hd, dp, dq
  integer(8) :: ntouch(LMAX), nsame(LMAX), ncount, n3touch, n3same
  real(8) :: dsum(LMAX)

  read(*, nml=run)
  N = L * L
  allocate(sp(N), nbr(4, N), stk(N))
  do iy = 0, L - 1
     do ix = 0, L - 1
        i = 1 + ix + L * iy
        nbr(1, i) = 1 + modulo(ix + 1, L) + L * iy        ! +x
        nbr(2, i) = 1 + modulo(ix - 1, L) + L * iy        ! -x
        nbr(3, i) = 1 + ix + L * modulo(iy + 1, L)        ! +y
        nbr(4, i) = 1 + ix + L * modulo(iy - 1, L)        ! -y
     end do
  end do
  padd = 1.0d0 - exp(-2.0d0 * K)
  w4 = exp(-4.0d0 * K); w8 = exp(-8.0d0 * K)

  ! ---------------- generator initialisation
  select case (src)
  case (1:9)
     call rng_init(src, seed)
  case (S_GFSR, S_DEC)
     call gfsr_init(g1, lags, ntap, seed, .false.)
  case (S_XOR2, S_ADD2)
     call gfsr_init(g1, lags, ntap, seed, .false.)
     call gfsr_init(g2, lags, ntap, seed2, .false.)
  case (S_ALFG)
     call gfsr_init(g1, lags, ntap, seed, .true.)
  case (S_GMT)
     call gfsr_init(g1, lags, ntap, seed, .false.)
     call rng_init(G_MT, seed2)
  end select

  if (alg == 'check') then
     call check_mode()
     stop
  end if
  if (alg == 'diag') then
     call diag_mode()
     stop
  end if

  ! ---------------- initial state
  if (ordered_start) then
     sp = 1
  else
     do i = 1, N
        sp(i) = merge(1, -1, u() < 0.5d0)
     end do
  end if

  call cpu_time(t0)
  do t = 1, ntherm
     call one_update(csize)
  end do

  nper = int(nmeas / nbatch)
  acf = 0; ring = 0; it = 0
  open(10, file=trim(outfile), status='replace')
  write(10, '(a)') '# ising.f90 output'
  write(10, '(a,a,a,a,a,i0,a,i0,a,i0)') '# alg=', trim(alg), ' order=', trim(order), ' L=', L, ' src=', src, ' dec=', dec
  write(10, '(a,4(1x,i0),a,i0,a,i0,a,i0)') '# lags=', lags, ' ntap=', ntap, ' seed=', seed, ' seed2=', seed2
  write(10, '(a,es24.16,a,i0,a,i0,a,i0,a,i0)') '# K=', K, ' ntherm=', ntherm, ' nbatch=', nbatch, ' nper=', nper, ' N=', N
  write(10, '(a)') '# batch  n  sumB  sumB2  sumCsize  sumAbsM  sumM2   (B = sum_<ij> s_i s_j, M = sum s_i)'
  do ib = 1, nbatch
     sB = 0; sB2 = 0; sC = 0; sAM = 0; sM2 = 0
     do t = 1, nper
        call one_update(csize)
        call measure(B, M)
        sB = sB + B; sB2 = sB2 + B * B; sC = sC + csize
        sAM = sAM + abs(M); sM2 = sM2 + M * M
        ! energy autocorrelation sums: acf(k) = sum_t B_t B_{t-k}
        ring(mod(it, int(NACF, 8))) = B
        if (it >= NACF) then
           do i = 0, NACF - 1
              acf(i) = acf(i) + B * ring(mod(it - i + NACF, int(NACF, 8)))
           end do
        end if
        it = it + 1
     end do
     write(10, '(i6,1x,i0,5(1x,i0))') ib, nper, sB, sB2, sC, sAM, sM2
  end do
  call cpu_time(t1)
  write(10, '(a,i0)') '# acf_nterms ', it - NACF
  do i = 0, NACF - 1
     write(10, '(a,i3,1x,i0)') '# acf ', i, acf(i)
  end do
  write(10, '(a,i0,a,f12.3)') '# ndraw ', ndraw, '  cpu_s ', t1 - t0
  close(10)

contains

  subroutine one_update(cs)
    integer, intent(out) :: cs
    select case (alg)
    case ('wolff')
       call wolff(cs)
    case ('metro')
       call metro_seq(); cs = 0
    case ('mrand')
       call metro_rand(); cs = 0
    case ('mskip')
       ! typewriter sweep plus one discarded draw: N+1 draws per sweep, so the
       ! same site sees random numbers N+1 (not N) draws apart
       call metro_seq(); cs = 0
       if (u() < -1.0d0) stop 'never'
    case default
       stop 'unknown alg'
    end select
  end subroutine one_update

  subroutine wolff(cs)
    integer, intent(out) :: cs
    integer :: i0, s0, j, kk, nb, top, head
    i0 = int(u() * N) + 1; if (i0 > N) i0 = N
    s0 = sp(i0); sp(i0) = -s0
    stk(1) = i0
    if (order == 'lifo') then
       top = 1; cs = 1
       do while (top > 0)
          j = stk(top); top = top - 1
          do kk = 1, 4
             nb = nbr(kk, j)
             if (sp(nb) == s0) then
                if (u() < padd) then
                   sp(nb) = -s0; top = top + 1; stk(top) = nb; cs = cs + 1
                end if
             end if
          end do
       end do
    else
       head = 1; top = 1
       do while (head <= top)
          j = stk(head); head = head + 1
          do kk = 1, 4
             nb = nbr(kk, j)
             if (sp(nb) == s0) then
                if (u() < padd) then
                   sp(nb) = -s0; top = top + 1; stk(top) = nb
                end if
             end if
          end do
       end do
       cs = top
    end if
  end subroutine wolff

  ! typewriter sweep, one random number drawn for every site
  subroutine metro_seq()
    integer :: i, h, de
    real(8) :: r
    do i = 1, N
       h = sp(nbr(1, i)) + sp(nbr(2, i)) + sp(nbr(3, i)) + sp(nbr(4, i))
       de = 2 * sp(i) * h          ! change of B is -de ; energy change is +de
       r = u()
       if (de <= 0) then
          sp(i) = -sp(i)
       else if (de == 4) then
          if (r < w4) sp(i) = -sp(i)
       else
          if (r < w8) sp(i) = -sp(i)
       end if
    end do
  end subroutine metro_seq

  ! N random-site attempts per sweep, two random numbers per attempt
  subroutine metro_rand()
    integer :: n1, i, h, de
    real(8) :: r
    do n1 = 1, N
       i = int(u() * N) + 1; if (i > N) i = N
       h = sp(nbr(1, i)) + sp(nbr(2, i)) + sp(nbr(3, i)) + sp(nbr(4, i))
       de = 2 * sp(i) * h
       r = u()
       if (de <= 0) then
          sp(i) = -sp(i)
       else if (de == 4) then
          if (r < w4) sp(i) = -sp(i)
       else
          if (r < w8) sp(i) = -sp(i)
       end if
    end do
  end subroutine metro_rand

  subroutine measure(Bo, Mo)
    integer(8), intent(out) :: Bo, Mo
    integer :: i, b, m
    b = 0; m = 0
    do i = 1, N
       b = b + sp(i) * (sp(nbr(1, i)) + sp(nbr(3, i)))
       m = m + sp(i)
    end do
    Bo = b; Mo = m
  end subroutine measure

  ! ---------------------------------------------------------------
  ! check: (a) own GFSR(103,250) reproduces rnglib R250 bit for bit;
  ! (b) three-point correlation T(p,q) = <(1/2-u_n)(1/2-u_{n-p})(1/2-u_{n-q})>
  !     of the selected source, with lags(1), lags(ntap), from nmeas draws.
  subroutine check_mode()
    integer(8), parameter :: NCHK = 10000000_8
    integer(8) :: n, nbad, ww
    real(8), allocatable :: hist(:)
    real(8) :: s3, s3sq, a, prod, kap, pacc_rr, pacc_aa, pacc_ar
    integer :: p, q, ia, ib1, ic
    integer(8) :: c111(0:1, 0:1, 0:1)
    type(gfsr_t) :: gt
    call rng_init(G_R250, seed)
    call gfsr_init(gt, [103, 250, 0, 0], 2, seed, .false.)
    nbad = 0
    do n = 1, NCHK
       ww = gfsr_word(gt)
       if (abs(dble(ww) / TWO32 - rng_next(G_R250)) > 0.0d0) nbad = nbad + 1
    end do
    write(*, '(a,i0,a,i0)') 'own GFSR(103,250) vs rnglib R250: mismatches = ', nbad, ' of ', NCHK
    p = lags(1); q = lags(ntap)
    allocate(hist(0:q))
    do n = 0, q
       hist(n) = u()
    end do
    s3 = 0; s3sq = 0; c111 = 0
    do n = q + 1, q + nmeas
       a = u()
       hist(mod(n, int(q + 1, 8))) = a
       prod = (0.5d0 - a) * (0.5d0 - hist(mod(n - p, int(q + 1, 8)))) * (0.5d0 - hist(mod(n - q, int(q + 1, 8))))
       s3 = s3 + prod; s3sq = s3sq + prod * prod
       ia = merge(1, 0, a < padd); ib1 = merge(1, 0, hist(mod(n - p, int(q + 1, 8))) < padd)
       ic = merge(1, 0, hist(mod(n - q, int(q + 1, 8))) < padd)
       c111(ia, ib1, ic) = c111(ia, ib1, ic) + 1
    end do
    ! connected third cumulant of the acceptance indicators I = [u < p_add]
    kap = 0
    do ia = 0, 1
       do ib1 = 0, 1
          do ic = 0, 1
             kap = kap + dble(c111(ia, ib1, ic)) * (ia - padd) * (ib1 - padd) * (ic - padd)
          end do
       end do
    end do
    kap = kap / dble(nmeas)
    pacc_rr = dble(c111(1, 0, 0)) / dble(c111(1, 0, 0) + c111(0, 0, 0))
    pacc_aa = dble(c111(1, 1, 1)) / dble(c111(1, 1, 1) + c111(0, 1, 1))
    pacc_ar = dble(c111(1, 1, 0) + c111(1, 0, 1)) / dble(c111(1, 1, 0) + c111(1, 0, 1) + c111(0, 1, 0) + c111(0, 0, 1))
    write(*, '(a,es13.5,a,es10.3,a,3f10.6)') '   kappa3(I) = ', kap, ' +- ', &
         sqrt((padd * (1 - padd))**3 / dble(nmeas)), '   P(acc | rej,rej), P(acc | acc,acc), P(acc | mixed) = ', &
         pacc_rr, pacc_aa, pacc_ar
    s3 = s3 / dble(nmeas); s3sq = s3sq / dble(nmeas)
    write(*, '(a,i0,a,i0,a,i0,a,i0,a,es14.6,a,es11.3,a,f12.2)') 'src=', src, ' dec=', dec, ' T(', p, ',', q, ') = ', s3, &
         ' +- ', sqrt((s3sq - s3 * s3) / dble(nmeas)), '   z = ', s3 / sqrt((s3sq - s3 * s3) / dble(nmeas))
  end subroutine check_mode

  ! ---------------------------------------------------------------
  ! diag: geometry of random-number consumption. For every draw store the
  ! target site and the cluster (or sweep) index; for lags l = 1..LMAX
  ! accumulate P(target sites within Manhattan distance <= 1) and
  ! P(same cluster). Also the joint probability for the pair (lags(1), lags(2)).
  subroutine diag_mode()
    integer :: lg, j, kk, nb, i0, s0, top, cl
    integer(8) :: nd0
    integer(8), allocatable :: hdraw(:)
    allocate(hdraw(0:4 * N + 1)); hdraw = 0
    dp = lags(1); dq = lags(2)
    ntouch = 0; nsame = 0; ncount = 0; n3touch = 0; n3same = 0; dsum = 0
    hd = 0; tsite = 1; tclus = -1
    do i0 = 1, N
       sp(i0) = merge(1, -1, u() < 0.5d0)
    end do
    do cl = 1, int(ntherm + nmeas)
       if (order == 'metr') then
          do j = 1, N
             call record(j, cl)
             if (u() < 0.0d0) stop
          end do
          cycle
       end if
       nd0 = ndraw
       i0 = int(u() * N) + 1; if (i0 > N) i0 = N
       call record(i0, cl)
       s0 = sp(i0); sp(i0) = -s0
       top = 1; stk(1) = i0
       do while (top > 0)
          j = stk(top); top = top - 1
          do kk = 1, 4
             nb = nbr(kk, j)
             if (sp(nb) == s0) then
                call record(nb, cl)
                if (u() < padd) then
                   sp(nb) = -s0; top = top + 1; stk(top) = nb
                end if
             end if
          end do
       end do
       if (cl > ntherm) hdraw(ndraw - nd0) = hdraw(ndraw - nd0) + 1
    end do
    open(11, file=trim(outfile), status='replace')
    write(11, '(a,i0,a,a,a,i0,a,i0)') '# diag L=', L, ' order=', order, ' nclusters=', nmeas, ' ndraws=', ncount
    write(11, '(a,i0,a,i0,a,es14.6,a,es14.6)') '# triple p=', dp, ' q=', dq, ' P3touch=', dble(n3touch)/dble(ncount), &
         ' P3same=', dble(n3same)/dble(ncount)
    write(11, '(a)') '# lag  P_touch  P_same  mean_dist'
    do lg = 1, LMAX
       write(11, '(i4,3(1x,es14.6))') lg, dble(ntouch(lg))/dble(ncount), dble(nsame(lg))/dble(ncount), dsum(lg)/dble(ncount)
    end do
    write(11, '(a)') '# draws-per-cluster histogram: d count'
    do lg = 0, 4 * N + 1
       if (hdraw(lg) > 0) write(11, '(a,i6,1x,i0)') '#H ', lg, hdraw(lg)
    end do
    close(11)
  end subroutine diag_mode

    subroutine record(site, c)
      integer, intent(in) :: site, c
      integer :: ll, hh, d, ip, iq
      hd = mod(hd + 1, LMAX + 1)
      tsite(hd) = site; tclus(hd) = c
      if (c <= ntherm) return
      ncount = ncount + 1
      do ll = 1, LMAX
         hh = mod(hd - ll + LMAX + 1, LMAX + 1)
         d = tdist(site, tsite(hh))
         if (d <= 1) ntouch(ll) = ntouch(ll) + 1
         if (tclus(hh) == c) nsame(ll) = nsame(ll) + 1
         dsum(ll) = dsum(ll) + d
      end do
      ip = mod(hd - dp + LMAX + 1, LMAX + 1); iq = mod(hd - dq + LMAX + 1, LMAX + 1)
      if (tdist(site, tsite(ip)) <= 1 .and. tdist(site, tsite(iq)) <= 1) n3touch = n3touch + 1
      if (tclus(ip) == c .and. tclus(iq) == c) n3same = n3same + 1
    end subroutine record
    integer function tdist(a, b)
      integer, intent(in) :: a, b
      integer :: dx, dy
      dx = abs(mod(a - 1, L) - mod(b - 1, L)); dx = min(dx, L - dx)
      dy = abs((a - 1) / L - (b - 1) / L);     dy = min(dy, L - dy)
      tdist = dx + dy
    end function tdist

end program ising
