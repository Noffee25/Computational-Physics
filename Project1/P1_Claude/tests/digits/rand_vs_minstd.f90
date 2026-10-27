! which truncation of the Park-Miller state reproduces gfortran RAND()?
!  A: floor(k * 2^22 / (2^31-1)) / 2^22      (MINSTD value truncated)
!  B: floor(k / 2^9) / 2^22 = (k with 9 low bits cleared) / 2^31
!  C: floor((k-1) / 2^9) / 2^22   (libgfortran: ((irand-1) & mask) / 2^31)
program rand_vs_minstd
  use rnglib
  implicit none
  integer(8) :: n, nd, k, ja, jb, jc, jr, misa, misb, misc
  real(8) :: xr, xm
  nd = 100000000_8
  call rng_init(G_RAND, 777_8); call rng_init(G_MINSTD, 777_8)
  misa = 0; misb = 0; misc = 0
  do n = 1, nd
     xr = rng_next(G_RAND); xm = rng_next(G_MINSTD); k = rng_state_int(G_MINSTD)
     jr = int(xr * 4194304.0d0, 8)
     ja = (k * 4194304_8) / 2147483647_8
     jb = k / 512_8
     jc = (k - 1_8) / 512_8
     if (jr /= jc) misc = misc + 1
     if (jr /= ja) misa = misa + 1
     if (jr /= jb) misb = misb + 1
  end do
  print *, 'N', nd, ' mismatches A (k/M truncated)', misa, ' B (k/2^31 truncated)', misb, ' C ((k-1)/2^31 truncated)', misc
end program rand_vs_minstd
