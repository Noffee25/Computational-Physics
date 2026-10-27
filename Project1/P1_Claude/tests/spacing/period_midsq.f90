! =====================================================================
!  period_midsq.f90 : complete functional-graph analysis of the
!  8-digit von Neumann middle-square map f(x) = floor(x^2/10^4) mod 10^8
!  on all 10^8 states (every possible seed).
!
!  For every state: tail length mu (steps until the orbit enters its
!  cycle) and the cycle it ends on. Memoized walk, each state visited
!  O(1) times. Writes
!     midsq_cycles.txt : cycle id, length, min element, basin size,
!                        mean and max tail of the basin
!     midsq_tail_hist.txt : histogram of mu over all seeds
!     midsq_rho_hist.txt  : histogram of mu+lambda (number of distinct
!                           values the generator ever produces)
!  Also verifies f against rnglib for 10^5 steps.
! =====================================================================
program period_midsq
  use rnglib
  implicit none
  integer, parameter :: NS = 100000000
  integer(2), allocatable :: tail(:)       ! -2 unvisited, -1 on current path, >=0 tail
  integer(2), allocatable :: cyc(:)        ! cycle id (1-based)
  integer, allocatable :: path(:)
  integer :: x0, x, np, i, k, pos, ncyc, c, maxtail
  integer, parameter :: MAXC = 30000, MAXP = 2000000
  integer(8) :: clen(MAXC), cmin(MAXC), basin(MAXC), tailsum(MAXC)
  integer :: tailmax(MAXC)
  integer(8), allocatable :: thist(:), rhist(:)
  integer(8) :: y, z
  integer(8) :: t0, t1, rate

  ! ---- verify the map against rnglib
  call rng_init(G_MIDSQ, 12345678_8)
  z = 12345678_8
  do i = 1, 100000
     y = int(rng_next(G_MIDSQ) * 0.0d0, 8)
     z = f(int(z))
     if (z /= rng_state_int(G_MIDSQ)) stop 'midsq map mismatch'
  end do
  print *, '# map verified'

  call system_clock(t0, rate)
  allocate(tail(0:NS-1), cyc(0:NS-1), path(MAXP))
  tail = -2_2
  cyc = 0_2
  ncyc = 0
  maxtail = 0
  do x0 = 0, NS - 1
     if (tail(x0) /= -2_2) cycle
     np = 0
     x = x0
     do while (tail(x) == -2_2)
        tail(x) = -1_2
        np = np + 1
        if (np > MAXP) stop 'path too long'
        path(np) = x
        x = int(f(x))
     end do
     if (tail(x) == -1_2) then
        ! new cycle: x is on the current path
        pos = 0
        do k = np, 1, -1
           if (path(k) == x) then
              pos = k; exit
           end if
        end do
        ncyc = ncyc + 1
        if (ncyc > MAXC) stop 'too many cycles'
        clen(ncyc) = np - pos + 1
        cmin(ncyc) = huge(1_8)
        do k = pos, np
           tail(path(k)) = 0_2
           cyc(path(k)) = int(ncyc, 2)
           cmin(ncyc) = min(cmin(ncyc), int(path(k), 8))
        end do
        do k = pos - 1, 1, -1
           if (pos - k > 32000) stop 'tail overflow'
           tail(path(k)) = int(pos - k, 2)
           cyc(path(k)) = int(ncyc, 2)
        end do
     else
        ! reached a state with known tail
        c = cyc(x)
        do k = np, 1, -1
           if (int(tail(x)) + (np - k + 1) > 32000) stop 'tail overflow'
           tail(path(k)) = int(int(tail(x)) + (np - k + 1), 2)
           cyc(path(k)) = int(c, 2)
        end do
     end if
  end do
  call system_clock(t1)
  print '(a,f8.2,a,i0)', '# graph done in ', dble(t1 - t0) / dble(rate), ' s; cycles = ', ncyc

  ! ---- statistics per cycle and histograms
  basin = 0; tailsum = 0; tailmax = 0
  maxtail = int(maxval(tail))
  allocate(thist(0:maxtail), rhist(0:maxtail + 70000))
  thist = 0; rhist = 0
  do x = 0, NS - 1
     c = cyc(x)
     basin(c) = basin(c) + 1
     tailsum(c) = tailsum(c) + tail(x)
     tailmax(c) = max(tailmax(c), int(tail(x)))
     thist(tail(x)) = thist(tail(x)) + 1
     k = int(tail(x)) + int(clen(c))
     if (k <= ubound(rhist, 1)) rhist(k) = rhist(k) + 1
  end do
  open(20, file='midsq_cycles.txt', status='replace')
  write(20, '(a)') '# id length min_element basin mean_tail max_tail'
  do c = 1, ncyc
     write(20, '(i0,1x,i0,1x,i0,1x,i0,1x,f12.4,1x,i0)') c, clen(c), cmin(c), basin(c), &
          dble(tailsum(c)) / dble(basin(c)), tailmax(c)
  end do
  close(20)
  open(21, file='midsq_tail_hist.txt', status='replace')
  write(21, '(a)') '# mu count'
  do k = 0, maxtail
     if (thist(k) > 0) write(21, '(i0,1x,i0)') k, thist(k)
  end do
  close(21)
  open(22, file='midsq_rho_hist.txt', status='replace')
  write(22, '(a)') '# mu_plus_lambda count'
  do k = 0, ubound(rhist, 1)
     if (rhist(k) > 0) write(22, '(i0,1x,i0)') k, rhist(k)
  end do
  close(22)
  ! a few named seeds
  open(23, file='midsq_seeds.txt', status='replace')
  write(23, '(a)') '# seed mu cycle_length cycle_min'
  do i = 1, 6
     select case (i)
     case (1); x = 777
     case (2); x = 12345678
     case (3); x = 65536
     case (4); x = 314159
     case (5); x = 77777777
     case (6); x = 99999999
     end select
     write(23, '(i0,1x,i0,1x,i0,1x,i0)') x, tail(x), clen(cyc(x)), cmin(cyc(x))
  end do
  close(23)

contains

  pure function f(x) result(y)
    integer, intent(in) :: x
    integer(8) :: y
    y = modulo((int(x, 8) * int(x, 8)) / 10000_8, 100000000_8)
  end function f

end program period_midsq
