! verifier: gfortran RAND() vs floor((s-1)/2^9)/2^22 over one FULL MINSTD period,
! with s computed here in integer(8) (no rnglib). Also checks the alternative rule.
program rand_full
  implicit none
  integer(8), parameter :: M = 2147483647_8, P = 2147483646_8
  integer(8) :: s, i, k, nbad, nbad2, k2, cnt_top, cnt_zero
  real(4) :: r
  call srand(777)
  s = 777_8; nbad = 0; nbad2 = 0; cnt_top = 0; cnt_zero = 0
  do i = 1, P
     r = rand()
     s = mod(16807_8 * s, M)
     k = (s - 1_8) / 512_8
     k2 = (s * 4194304_8) / M
     if (int(dble(r) * 4194304.0d0, 8) /= k .or. dble(r) * 4194304.0d0 /= dble(k)) nbad = nbad + 1
     if (dble(r) * 4194304.0d0 /= dble(k2)) nbad2 = nbad2 + 1
     if (k == 4194303_8) cnt_top = cnt_top + 1
     if (k == 0_8) cnt_zero = cnt_zero + 1
  end do
  print '(A,I12,A,I10,A,I10,A,I6,A,I6,A,I12)', 'draws', P, ' mismatch(s-1)/512:', nbad, &
       ' mismatch 2^22 s/M:', nbad2, ' n(top)=', cnt_top, ' n(0)=', cnt_zero, ' final state', s
end program rand_full
