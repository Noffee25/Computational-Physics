! =====================================================================
!  lagscan.f90 : blind search for three-term XOR identities
!
!  usage: lagscan gid seed T L outfile
!
!  w(n) = floor(x_n 2^32) (top 32 bits of each output).  For every lag
!  pair 1 <= s < t <= T the statistic
!       F(s,t) = sum_{n} popcount( w(n) XOR w(n-s) XOR w(n-t) )
!  over L consecutive n counts how often the XOR of three bits is 1.
!  Under the null (independent fair bits) F ~ Binomial(32 L, 1/2), so
!       z = (F - 16 L) / sqrt(8 L)  ~  N(0,1).
!  Bits that are constant over the whole sample (dead or stuck bits, e.g.
!  bits 23..32 of RAND, bits 29 and 31 of RANDU) are masked out, and the
!  null becomes Binomial(nb L, 1/2) with nb the number of live bits.
!  The program writes every pair with |z| > 5, the histogram of all z,
!  and the extreme values.
! =====================================================================
program lagscan
  use rnglib
  implicit none
  integer :: gid, tmax, s, t, ib, nbin
  integer(8) :: seed, lw, n, f, npair, hist(-41:41), wand, wor, mask
  integer :: nb
  integer(8), allocatable :: w(:), u(:)
  real(8) :: z, zmin, zmax, x
  character(len=256) :: arg, outf
  call get_command_argument(1, arg); read(arg, *) gid
  call get_command_argument(2, arg); read(arg, *) seed
  call get_command_argument(3, arg); read(arg, *) tmax
  call get_command_argument(4, arg); read(arg, *) lw
  call get_command_argument(5, outf)
  call rng_init(gid, seed)
  allocate(w(0:lw + tmax - 1), u(0:lw - 1))
  do n = 0, lw + tmax - 1
     x = rng_next(gid)
     w(n) = int(x * 4294967296.0d0, 8)
  end do
  wand = 4294967295_8; wor = 0
  do n = 0, lw + tmax - 1
     wand = iand(wand, w(n)); wor = ior(wor, w(n))
  end do
  mask = iand(4294967295_8, ieor(wand, wor))      ! bits that change at least once
  nb = popcnt(mask)
  open(10, file=trim(outf), status='replace')
  write(10, '(a,a,a,i12,a,i6,a,i10,a,i3)') '# ', trim(gen_name(gid)), ' seed ', seed, ' T ', tmax, ' L ', lw, &
       ' livebits ', nb
  hist = 0; npair = 0; zmin = 1d300; zmax = -1d300
  do s = 1, tmax - 1
     do n = 0, lw - 1
        u(n) = ieor(w(n + tmax), w(n + tmax - s))
     end do
     do t = s + 1, tmax
        f = 0
        do n = 0, lw - 1
           f = f + popcnt(iand(ieor(u(n), w(n + tmax - t)), mask))
        end do
        z = (dble(f) - 0.5d0 * nb * dble(lw)) / sqrt(0.25d0 * nb * dble(lw))
        npair = npair + 1
        zmin = min(zmin, z); zmax = max(zmax, z)
        ib = nint(z * 4.0d0)
        if (ib < -40) ib = -41
        if (ib > 40) ib = 41
        hist(ib) = hist(ib) + 1
        if (abs(z) > 5.0d0) write(10, '(a,2i6,i12,f14.3)') 'HIT ', s, t, f, z
     end do
  end do
  write(10, '(a,i12,2f14.3)') 'NPAIR ', npair, zmin, zmax
  do ib = -41, 41
     write(10, '(a,f8.3,i12)') 'HIST ', ib / 4.0d0, hist(ib)
  end do
  close(10)
  write(*, '(a,a,a,i10,2f12.3)') 'done ', trim(gen_name(gid)), ' pairs ', npair, zmin, zmax
end program lagscan
