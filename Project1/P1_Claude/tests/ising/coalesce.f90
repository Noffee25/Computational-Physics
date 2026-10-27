! Coalescence of two Metropolis (typewriter) chains driven by common random
! numbers: the mechanism that makes a chain driven by a generator of period P
! (P a multiple of the sweep length N) exactly periodic once X_t = X_{t+P}.
! Chains start from independent random states (all-up versus all-down would
! be mirror images forever by the global spin-flip symmetry); report the
! number of sweeps until identical (kind=1) or mirror images (kind=-1).
program coalesce
  use rnglib
  implicit none
  integer, parameter :: L = 16, N = L * L, NTRIAL = 20
  integer(8), parameter :: TMAX = 200000_8
  integer :: sa(N), sb(N), nbr(4, N), i, ix, iy, trial, ha, hb, dea, deb
  integer(8) :: t
  integer :: kind
  real(8) :: K, w4, w8, r
  K = 0.5d0 * log(1.0d0 + sqrt(2.0d0)); w4 = exp(-4 * K); w8 = exp(-8 * K)
  do iy = 0, L - 1
     do ix = 0, L - 1
        i = 1 + ix + L * iy
        nbr(1, i) = 1 + modulo(ix + 1, L) + L * iy; nbr(2, i) = 1 + modulo(ix - 1, L) + L * iy
        nbr(3, i) = 1 + ix + L * modulo(iy + 1, L); nbr(4, i) = 1 + ix + L * modulo(iy - 1, L)
     end do
  end do
  call rng_init(G_MT, 777_8)
  do trial = 1, NTRIAL
     kind = 0
     do i = 1, N
        sa(i) = merge(1, -1, rng_next(G_MT) < 0.5d0); sb(i) = merge(1, -1, rng_next(G_MT) < 0.5d0)
     end do
     do t = 1, TMAX
        do i = 1, N
           r = rng_next(G_MT)
           ha = sa(nbr(1,i)) + sa(nbr(2,i)) + sa(nbr(3,i)) + sa(nbr(4,i)); dea = 2 * sa(i) * ha
           hb = sb(nbr(1,i)) + sb(nbr(2,i)) + sb(nbr(3,i)) + sb(nbr(4,i)); deb = 2 * sb(i) * hb
           if (dea <= 0) then; sa(i) = -sa(i); else if (dea == 4) then; if (r < w4) sa(i) = -sa(i); else; if (r < w8) sa(i) = -sa(i); end if
           if (deb <= 0) then; sb(i) = -sb(i); else if (deb == 4) then; if (r < w4) sb(i) = -sb(i); else; if (r < w8) sb(i) = -sb(i); end if
        end do
        if (all(sa == sb)) then
           kind = 1; exit
        end if
        if (all(sa == -sb)) then
           kind = -1; exit
        end if
     end do
     write(*, '(i4,1x,i0,1x,i0)') trial, t, kind
  end do
end program coalesce
