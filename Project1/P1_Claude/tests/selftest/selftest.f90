program selftest
  use rnglib
  implicit none
  integer :: i, id, k
  integer(8) :: w(0:300), y
  real(8) :: x, xs(3), t0, t1, s
  logical :: ok
  ! 1. RANDU seed 1
  call rng_init(G_RANDU, 1_8)
  do i = 1, 3
     x = rng_next(G_RANDU); print '(a,i12)', 'RANDU state ', rng_state_int(G_RANDU)
  end do
  ! 2. RANDU 314159 triplets
  call rng_init(G_RANDU, 314159_8)
  do k = 1, 3
     do i = 1, 3; xs(i) = rng_next(G_RANDU); end do
     print '(a,3f9.5,f12.6)', 'RANDU trip ', xs, 9*xs(1)-6*xs(2)+xs(3)
  end do
  ! 3. MINSTD
  call rng_init(G_MINSTD, 1_8)
  do i = 1, 10000; x = rng_next(G_MINSTD); end do
  print '(a,i12,a)', 'MINSTD z10000 ', rng_state_int(G_MINSTD), ' (expect 1043618065)'
  ! 4. MT
  call rng_init(G_MT, 5489_8)
  y = mt_u32(); print '(a,i12,a)', 'MT first u32 ', y, ' (expect 3499211612)'
  y = mt_u32(); print '(a,i12,a)', 'MT second u32 ', y, ' (expect 581869302)'
  ! 5. R250 recurrence
  call rng_init(G_R250, 777_8)
  do i = 0, 300; w(i) = nint(rng_next(G_R250) * 4294967296.0d0, 8); end do
  ok = .true.
  do i = 250, 300
     if (ieor(ieor(w(i), w(i-103)), w(i-250)) /= 0) ok = .false.
  end do
  print *, 'R250 xor recurrence holds: ', ok
  ! 6. RAND vs MINSTD
  call rng_init(G_RAND, 777_8)
  call rng_init(G_MINSTD, 777_8)
  do i = 1, 5
     x = rng_next(G_RAND); s = rng_next(G_MINSTD)
     print '(a,f14.10,a,f14.10,a,f16.3,a,i12)', 'RAND ', x, ' MINSTD ', s, ' RAND*2^22 ', x*4194304d0, ' pm ', rng_state_int(G_MINSTD)
  end do
  ! 7. speed and mean
  do id = 1, NGEN
     call rng_init(id, 777_8)
     call cpu_time(t0); s = 0
     do i = 1, 100000000; s = s + rng_next(id); end do
     call cpu_time(t1)
     print '(a8,a,f8.5,a,f6.2,a)', gen_name(id), ' mean ', s/1d8, '  ', (t1-t0), ' s per 1e8'
  end do
end program
