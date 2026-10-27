! =====================================================================
! float_weights.f90 : bin probabilities after rounding a continuous
! uniform variate to REAL*4, exactly as rnglib does with rng_single.
!
! A REAL*4 value f receives every double x in its round-to-nearest cell
! [(f_prev+f)/2, (f+f_next)/2); the largest float below 1 also receives
! [1-2^-25, 1) because rnglib maps the rounded value 1.0 back to it.
! The bin probability p_b is the sum of the cell widths of the floats
! whose value lies in bin b = int(f*B)+1.  Floats are visited in
! increasing order, so each bin is completed before the next starts and
! no per-bin array is needed (streaming over ~3.4e8 floats).
!
! Output per depth (B = 10^d, d = 1..8): sum of delta_b^2 with
! delta_b = B p_b - 1, the number of bins with p_b = 0, and a histogram
! of w_b = B p_b in cells of width 1e-6 (nonzero cells only).
! usage: float_weights OUTPREFIX
! =====================================================================
program float_weights
  implicit none
  integer, parameter :: NW = 8000000
  integer(8), allocatable :: wh(:)
  integer :: d, nb, bcur, b, ibits, ibot, itop, nzero, j, k
  real(8) :: db, fprev, f, fnext, lo, hi, w, acc, sdel2, sumw
  real(4) :: r4
  character(len=128) :: pre
  call get_command_argument(1, pre)
  allocate(wh(0:NW-1))
  open(10, file=trim(pre)//'_w.out', status='unknown')
  open(11, file=trim(pre)//'_sum.out', status='unknown')
  ibot = transfer(2.0**(-40), ibits)      ! bit pattern of 2^-40
  itop = transfer(1.0, ibits) - 1         ! largest float below 1
  do d = 1, 8
     nb = 10**d; db = dble(nb)
     wh = 0; sdel2 = 0.0d0; nzero = 0; sumw = 0.0d0
     bcur = 1; acc = 0.0d0
     fprev = 0.0d0
     do ibits = ibot, itop
        f = dble(transfer(ibits, r4))
        if (ibits < itop) then
           fnext = dble(transfer(ibits + 1, r4))
           hi = 0.5d0 * (f + fnext)
        else
           hi = 1.0d0
        end if
        if (ibits == ibot) then
           lo = 0.0d0               ! continuous mass below 2^-40 joins bin 1
        else
           lo = 0.5d0 * (fprev + f)
        end if
        w = hi - lo
        b = int(f * db) + 1
        if (b > nb) b = nb
        if (b /= bcur) then
           call flush_bin(acc)
           do j = bcur + 1, b - 1        ! bins holding no float at all
              call flush_bin(0.0d0)
           end do
           bcur = b; acc = 0.0d0
        end if
        acc = acc + w
        fprev = f
     end do
     call flush_bin(acc)
     do j = bcur + 1, nb
        call flush_bin(0.0d0)
     end do
     write(11, '(I1,1X,I10,1X,ES22.14,1X,I10,1X,ES22.14)') d, nb, sdel2, nzero, sumw
     do k = 0, NW - 1
        if (wh(k) > 0) write(10, '(I1,1X,F12.7,1X,I11)') d, (dble(k) + 0.5d0) * 1.0d-6, wh(k)
     end do
  end do
  close(10); close(11)
contains
  subroutine flush_bin(pw)
    real(8), intent(in) :: pw
    real(8) :: wb
    integer :: kk
    wb = pw * db
    sumw = sumw + pw
    sdel2 = sdel2 + (wb - 1.0d0)**2
    if (pw == 0.0d0) nzero = nzero + 1
    kk = int(wb * 1.0d6)
    if (kk > NW - 1) kk = NW - 1
    wh(kk) = wh(kk) + 1
  end subroutine flush_bin
end program float_weights
