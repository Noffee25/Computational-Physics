! exact full-period fraction of RANDU triplets (mapped to [-1,1]^3) inside the unit
! ball that lie on the central plane 9x-6y+z = 0 (i.e. 9X-6Y+Z = 2^32 in integers)
program lat_ballfrac
  implicit none
  integer(8), parameter :: M = 2147483648_8
  integer(8) :: x1, x2, x3, n, nin, ncen, nband
  real(8) :: v(3), r2, ct
  x1 = 777; x2 = iand(65539_8 * x1, M - 1)
  nin = 0; ncen = 0; nband = 0
  do n = 1, 536870912_8
     x3 = iand(65539_8 * x2, M - 1)
     v = 2d0 * [dble(x1), dble(x2), dble(x3)] / dble(M) - 1d0
     r2 = sum(v * v)
     if (r2 <= 1d0 .and. r2 > 0d0) then
        nin = nin + 1
        if (9 * x1 - 6 * x2 + x3 == 2 * M) ncen = ncen + 1
        ct = (9d0 * v(1) - 6d0 * v(2) + v(3)) / sqrt(118d0 * r2)
        if (abs(ct) > 1d-9 .and. abs(ct) < 2d0 / sqrt(118d0) - 1d-9) nband = nband + 1
     end if
     x1 = x2; x2 = x3
  end do
  print '(a,i12,a,i12,a,f12.8,a,i6)', 'in ball ', nin, '  on central plane ', ncen, '  fraction ', &
       dble(ncen) / dble(nin), '  in forbidden band ', nband
end program
