! =====================================================================
!  cells.f90 : d-dimensional serial (cell occupancy) test
!  Non-overlapping d-tuples are dropped into k^d equal cells for several
!  resolutions k at once.  For each k we write
!      sum_i n_i^2 (exact), number of empty cells, max count,
!  from which X^2 = sum n_i^2/lambda - Nt  (lambda = Nt/k^d).
!  mode 0 : tuples (x_{t}, x_{t+1}, ..., x_{t+d-1}) from consecutive values
!  mode 1 : lagged triples (x_n, x_{n-103}, x_{n-250}) built from disjoint
!           elements: in each window of 353 values, j = 0..102 gives the
!           triple (w(j), w(j+147), w(j+250)), which uses every element at
!           most once (gaps 103 and 147 are both >= 103).
!  usage: cells gen seed d mode Nt nk k1 k2 ... outfile
! =====================================================================
program cells
  use rnglib
  implicit none
  integer, parameter :: MAXK = 16
  integer :: gen, d, mode, nk, kk(MAXK), ik, i, jj, c
  integer(8) :: seed, Nt, t, idx, ncell(MAXK), off(MAXK), tot, sq(MAXK), nempty(MAXK), mx(MAXK)
  integer, allocatable :: h(:)
  real(8) :: x(8), w(0:352), t0, t1
  character(len=256) :: arg, outfile
  call get_command_argument(1, arg); read(arg,*) gen
  call get_command_argument(2, arg); read(arg,*) seed
  call get_command_argument(3, arg); read(arg,*) d
  call get_command_argument(4, arg); read(arg,*) mode
  call get_command_argument(5, arg); read(arg,*) Nt
  call get_command_argument(6, arg); read(arg,*) nk
  do ik = 1, nk
     call get_command_argument(6+ik, arg); read(arg,*) kk(ik)
  end do
  call get_command_argument(7+nk, outfile)
  tot = 0
  do ik = 1, nk
     ncell(ik) = int(kk(ik),8)**d
     off(ik) = tot
     tot = tot + ncell(ik)
  end do
  allocate(h(0:tot-1)); h = 0
  call rng_init(gen, seed)
  call cpu_time(t0)
  if (mode == 0) then
     do t = 1, Nt
        do i = 1, d; x(i) = rng_next(gen); end do
        call drop()
     end do
  else
     t = 0
     do while (t < Nt)
        do i = 0, 352; w(i) = rng_next(gen); end do
        do jj = 0, 102
           if (t >= Nt) exit
           x(1) = w(jj+250); x(2) = w(jj+147); x(3) = w(jj)
           call drop(); t = t + 1
        end do
     end do
  end if
  call cpu_time(t1)
  open(11, file=trim(outfile), status='replace')
  write(11,'(a,a,a,i12,a,i2,a,i2,a,i12,a,f8.2)') '# ', gen_name(gen), ' seed ', seed, ' d ', d, &
       ' mode ', mode, ' Nt ', Nt, ' cpu ', t1-t0
  write(11,'(a)') '# k  ncells  sum_n2  nempty  maxcount'
  do ik = 1, nk
     sq(ik) = 0; nempty(ik) = 0; mx(ik) = 0
     do idx = off(ik), off(ik) + ncell(ik) - 1
        c = h(idx)
        sq(ik) = sq(ik) + int(c,8)*c
        if (c == 0) nempty(ik) = nempty(ik) + 1
        if (c > mx(ik)) mx(ik) = c
     end do
     write(11,'(i6,i12,i20,i12,i12)') kk(ik), ncell(ik), sq(ik), nempty(ik), mx(ik)
  end do
  close(11)
contains
  subroutine drop()
    integer :: q, m
    integer(8) :: id
    do q = 1, nk
       id = 0
       do m = 1, d
          id = id * kk(q) + min(int(x(m) * kk(q)), kk(q) - 1)
       end do
       h(off(q) + id) = h(off(q) + id) + 1
    end do
  end subroutine
end program
