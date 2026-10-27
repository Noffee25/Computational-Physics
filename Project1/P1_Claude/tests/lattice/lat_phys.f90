! =====================================================================
! lat_phys.f90 : two physical consequences of the RANDU planes.
!
! (1) "X-ray diffraction of an ideal gas".  N = 1e5 particles are placed in
!     the periodic unit box at r_j = (x_{3j+1}, x_{3j+2}, x_{3j+3}).  The static
!     structure factor S(n) = |sum_j exp(2 pi i n.r_j)|^2 / N is computed for
!     every integer wave vector n in [-10,10]^2 x [0,10].  For an ideal gas
!     S(n) ~ Exp(1) for n /= 0; a lattice gives Bragg peaks S = N.
!
! (2) "Isotropic emission".  Directions are drawn by the standard cube
!     rejection method (v = 2u - 1, accept 0 < |v| <= 1, d = v/|v|) and,
!     as a control, by the 2D method (cos theta = 2u_1 - 1, phi = 2 pi u_2).
!     Histograms of cos(theta) relative to the RANDU plane normal
!     n = (9,-6,1)/sqrt(118) and an equal-area sky map in (phi, cos theta).
! =====================================================================
program lat_phys
  use rnglib
  implicit none
  integer, parameter :: NPART = 100000, KM = 10
  integer, parameter :: NDIR = 2000000, NBC = 400, NPH = 144, NCT = 72
  integer :: ig, ids(5), j, k1, k2, k3, unit, im, nmeth
  character(len=8) :: nm(5)
  complex(8) :: ex(-KM:KM), ey(-KM:KM), ez(0:KM), A(-KM:KM, -KM:KM, 0:KM), exy
  real(8) :: r(3), S(-KM:KM, -KM:KM, 0:KM, 5), t0, t1
  real(8), parameter :: TWOPI = 6.283185307179586d0
  real(8) :: v(3), rr, d(3), e1(3), e2(3), e3(3), ct, phi, u1, u2
  integer(8) :: hc(NBC, 4), sky(NPH, NCT, 4), nexact(4), ndraw(4)
  integer :: ib, ip, ic, im2

  ids = [G_RANDU, G_RANDM, G_MINSTD, G_MT, G_RANNUM]
  nm = [character(len=8) :: 'RANDU', 'RANDM', 'MINSTD', 'MT19937', 'RANNUM']

  ! ------------------------------------------------ (1) structure factor
  do ig = 1, 5
     call cpu_time(t0)
     call rng_init(ids(ig), 777_8)
     A = (0d0, 0d0)
     do j = 1, NPART
        r(1) = rng_next(ids(ig)); r(2) = rng_next(ids(ig)); r(3) = rng_next(ids(ig))
        do k1 = 0, KM
           ex(k1) = exp(cmplx(0d0, TWOPI * k1 * r(1), 8))
           ey(k1) = exp(cmplx(0d0, TWOPI * k1 * r(2), 8))
           ez(k1) = exp(cmplx(0d0, TWOPI * k1 * r(3), 8))
           ex(-k1) = conjg(ex(k1)); ey(-k1) = conjg(ey(k1))
        end do
        do k3 = 0, KM
           do k2 = -KM, KM
              exy = ey(k2) * ez(k3)
              do k1 = -KM, KM
                 A(k1, k2, k3) = A(k1, k2, k3) + ex(k1) * exy
              end do
           end do
        end do
     end do
     S(:, :, :, ig) = (dble(A)**2 + aimag(A)**2) / dble(NPART)
     call cpu_time(t1)
     print '(a8,a,es12.4,a,3i4,a,f6.1,a)', nm(ig), ' max S(n /= 0) = ', maxval(S(:, :, :, ig), &
          mask=(S(:, :, :, ig) < 0.999d0 * NPART)), '   S(9,-6,1) = ', 9, -6, 1, '  cpu ', t1 - t0, ' s'
     print '(a,es14.6)', '      S(9,-6,1) = ', S(9, -6, 1, ig)
  end do
  open(newunit=unit, file='structure_factor.txt', status='replace')
  write(unit, '(a)') '# n1 n2 n3  S_RANDU S_RANDM S_MINSTD S_MT19937 S_RANNUM   (N = 1e5 particles, seed 777)'
  do k3 = 0, KM
     do k2 = -KM, KM
        do k1 = -KM, KM
           write(unit, '(3i4,5es14.6)') k1, k2, k3, S(k1, k2, k3, :)
        end do
     end do
  end do
  close(unit)

  ! ------------------------------------------------ (2) isotropy
  e3 = [9d0, -6d0, 1d0] / sqrt(118d0)
  e1 = [2d0, 3d0, 0d0] / sqrt(13d0)
  e2 = [e3(2) * e1(3) - e3(3) * e1(2), e3(3) * e1(1) - e3(1) * e1(3), e3(1) * e1(2) - e3(2) * e1(1)]
  hc = 0; sky = 0; nexact = 0; ndraw = 0
  ! methods: 1 RANDU cube rejection, 2 MT cube rejection, 3 RANDU 2D method, 4 RANDM cube rejection
  do im = 1, 4
     select case (im)
     case (1, 3); ig = G_RANDU
     case (2); ig = G_MT
     case (4); ig = G_RANDM
     end select
     call rng_init(ig, 777_8)
     j = 0
     do while (j < NDIR)
        if (im == 3) then
           u1 = rng_next(ig); u2 = rng_next(ig); ndraw(im) = ndraw(im) + 2
           ct = 2d0 * u1 - 1d0; phi = TWOPI * u2
           rr = sqrt(max(0d0, 1d0 - ct * ct))
           d = [rr * cos(phi), rr * sin(phi), ct]
        else
           v(1) = 2d0 * rng_next(ig) - 1d0
           v(2) = 2d0 * rng_next(ig) - 1d0
           v(3) = 2d0 * rng_next(ig) - 1d0
           ndraw(im) = ndraw(im) + 3
           rr = sqrt(sum(v * v))
           if (rr > 1d0 .or. rr == 0d0) cycle
           d = v / rr
        end if
        j = j + 1
        ct = sum(d * e3)
        if (abs(ct) < 1d-12) nexact(im) = nexact(im) + 1
        ib = min(NBC, int((ct + 1d0) / 2d0 * NBC) + 1)
        hc(ib, im) = hc(ib, im) + 1
        phi = atan2(sum(d * e2), sum(d * e1))
        ip = min(NPH, int((phi + 0.5d0 * TWOPI) / TWOPI * NPH) + 1)
        ic = min(NCT, int((ct + 1d0) / 2d0 * NCT) + 1)
        sky(ip, ic, im) = sky(ip, ic, im) + 1
     end do
     print '(a,i2,a,i9,a,f9.6,a,i12)', 'isotropy method', im, ': directions with |cos theta| < 1e-12: ', &
          nexact(im), '  fraction ', dble(nexact(im)) / NDIR, '  uniforms used ', ndraw(im)
  end do
  open(newunit=unit, file='isotropy_costheta.txt', status='replace')
  write(unit, '(a)') '# bin_center  RANDU_cube  MT_cube  RANDU_2D  RANDM_cube   (2e6 directions, 400 bins)'
  do ib = 1, NBC
     write(unit, '(f10.5,4i10)') -1d0 + (ib - 0.5d0) * 2d0 / NBC, hc(ib, :)
  end do
  close(unit)
  open(newunit=unit, file='isotropy_sky.txt', status='replace')
  write(unit, '(a)') '# ip ic  counts for 4 methods  (phi bins 144, cos theta bins 72)'
  do ic = 1, NCT
     do ip = 1, NPH
        write(unit, '(2i4,4i8)') ip, ic, sky(ip, ic, :)
     end do
  end do
  close(unit)
  open(newunit=unit, file='isotropy_exact.txt', status='replace')
  write(unit, '(a)') '# method  n_exact_on_great_circle  NDIR'
  do im = 1, 4
     write(unit, '(i3,2i10)') im, nexact(im), NDIR
  end do
  close(unit)
end program lat_phys
