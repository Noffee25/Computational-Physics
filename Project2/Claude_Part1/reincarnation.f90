!=======================================================================
!  reincarnation.f90  --  Conway's Game of Life, over and over again
!
!  Random "soups" are poured into a walled Petri dish (cells beyond the
!  rim are permanently dead) and every game is run until it STAGNATES,
!  i.e. the dish repeats an earlier state: it has died out, frozen into
!  still lifes or locked into oscillators.  The generation at which that
!  final cycle is first entered is the game's lifetime T.  Re-seeding the
!  random number generator for every game, we keep tallies of
!
!     1. T, extinction and left-over "ash" vs. the initial filling factor
!     2. T vs. the size of the Petri dish, at a fixed filling factor
!     3. the probability distribution of T
!     4. how big a seed colony must be to last 100+ generations
!
!  plus a census of how the games end.  While the batch runs, one of its
!  games is animated beside the tallies (cells are coloured by age, so
!  new growth glows and frozen debris fades); at the end the longest-lived
!  soup -- "the Methuselah" -- is replayed.
!
!  build:  gfortran -O2 reincarnation.f90 -o reincarnation
!  run:    ./reincarnation [seed] [fast]
!          needs a UTF-8, 256-colour terminal of at least 100 x 26.  "fast"
!          skips the animation and only prints the results, as does sending
!          the output to a file ( ./reincarnation > results.txt ).
!          If Ctrl-C leaves the colours odd, type  reset .
!=======================================================================
program reincarnation
   use, intrinsic :: iso_fortran_env, only: int8, int64, output_unit
   use, intrinsic :: iso_c_binding,   only: c_int
   implicit none

   interface        ! POSIX usleep() for smooth frames (SLEEP() only does whole
                    ! seconds) and isatty() to notice when output is redirected
      function usleep(usec) bind(c, name='usleep') result(ierr)
         import :: c_int
         integer(c_int), value :: usec
         integer(c_int)        :: ierr
      end function usleep
      function c_isatty(fd) bind(c, name='isatty') result(yes)
         import :: c_int
         integer(c_int), value :: fd
         integer(c_int)        :: yes
      end function c_isatty
   end interface

   ! ---------------------------------------------------------------- knobs
   integer, parameter :: LMAX  = 128          ! biggest dish
   integer, parameter :: PMAX  = 120          ! longest period looked for
   integer, parameter :: GMAX  = 50000        ! still going after this: "unsettled"
   integer, parameter :: TLONG = 100          ! long-lived, for experiment 4
   integer, parameter :: NT3   = 2000         ! games in the lifetime distribution
   integer, parameter :: SIZES(9) = [8, 12, 16, 24, 32, 48, 64, 96, 128]
   integer, parameter :: CW    = 40           ! the dish canvas shows CW x CW cells
   integer, parameter :: PW    = 54           ! inner width of the tally panel
   integer, parameter :: NROW  = CW/2 + 2     ! height of both boxes, in text rows
   integer, parameter :: TW    = CW + PW + 6  ! width of the whole dashboard
   integer, parameter :: LW    = 2048         ! storage length of one screen line
   real(8), parameter :: FPS   = 30d0         ! frames (= shown generations) per second

   ! -------------------------------------------------------------- looks
   character(*), parameter :: ESC = achar(27), NL = achar(10)
   character(len=:), allocatable :: RST, BOLD   ! left empty when output is redirected
   integer, parameter :: C_AGAR = 235, C_RIM = 240, C_DIM = 245, C_TITLE = 117, &
                         C_BAR = 79, C_HOT = 215
   integer, parameter :: C_AGE(6) = [230, 159, 116, 73, 37, 30]   ! newborn ... old
   character(len=3), parameter :: H8(8) = [character(len=3) :: &
                         '▏', '▎', '▍', '▌', '▋', '▊', '▉', '█']   ! eighth blocks
   character(len=3), parameter :: V8(8) = [character(len=3) :: &
                         '▁', '▂', '▃', '▄', '▅', '▆', '▇', '█']

   ! -------------------------------------------------------------- types
   type cond_t                                ! one experimental condition
      integer :: L = 40, side = 0, trials = 0 ! side > 0: side x side seed colony
      real(8) :: fill = 0.35d0
   end type cond_t

   type tally_t                               ! what is counted per condition
      integer :: n = 0, extinct = 0, osc = 0, long = 0, stuck = 0
      integer :: tmax = -1, seedmax = 0
      real(8) :: st = 0, st2 = 0, ash = 0
   end type tally_t

   type exp_t
      character(len=64) :: name = '', note = ''
      integer :: nc = 0
      type(cond_t)  :: c(19)
      type(tally_t) :: t(19)
   end type exp_t

   type game_t
      integer :: L = 0, side = 0, seed = 0, gen = 0, life = 0, period = 0
      real(8) :: fill = 0
      logical :: done = .false.
      integer(int8), allocatable :: c(:,:)    ! cells (0:L+1,0:L+1); the rim stays dead
      integer,       allocatable :: age(:,:)  ! generations each cell has been alive
      integer,       allocatable :: pop(:)    ! population history, pop(0:GMAX)
      integer(int64) :: h(0:PMAX-1) = 0       ! hashes of the last PMAX generations
   end type game_t

   ! ------------------------------------------------------------ globals
   type(exp_t)    :: ex(4)
   integer(int64) :: zkey(LMAX, LMAX)          ! Zobrist key of every cell
   integer        :: census(0:PMAX+1) = 0      ! 0 dead, 1 still, p period, PMAX+1 unsettled
   integer        :: life3(NT3) = 0, n3 = 0    ! lifetimes of experiment 3, kept sorted
   integer        :: base = 1000, nplayed = 0  ! base seed, games played so far
   integer        :: ngames = 0                ! games in a whole run
   integer        :: mspeed = 1                ! generations per frame in the replay
   logical        :: animate = .true., color = .true.
   real(8)        :: tstart
   character(len=262144) :: fb                 ! a whole frame, written in one go
   integer        :: fbn, ie

   call setup()
   tstart = now()
   if (animate) call emit(ESC//'[2J')
   do ie = 1, 4
      call run_experiment(ie)
   end do
   if (animate) call methuselah()
   call report()

contains

   ! ===================================================== setting up

   subroutine setup()
      character(len=32) :: arg
      integer(int64) :: j8
      integer :: i, j, ios
      real(8) :: u(2)
      do i = 1, command_argument_count()
         call get_command_argument(i, arg)
         if (arg == 'fast') then
            animate = .false.
            cycle
         end if
         ios = 1;  j8 = 0
         if (verify(trim(arg), '+-0123456789') == 0) read(arg, *, iostat=ios) j8
         if (ios /= 0 .or. abs(j8) > huge(base)) then
            write(output_unit, '(a)') 'usage: reincarnation [seed] [fast]      (seed: an integer)'
            stop
         end if
         base = int(j8)
      end do
      if (c_isatty(1_c_int) == 0) then         ! redirected: plain results, no animation
         animate = .false.
         color   = .false.
      end if
      RST = '';  BOLD = ''
      if (color) then
         RST = ESC//'[0m';  BOLD = ESC//'[1m'
      end if

      call seed_rng(271828)                    ! fixed keys for the state hash
      do j = 1, LMAX
         do i = 1, LMAX
            call random_number(u)
            zkey(i, j) = ior(ishft(int(u(1)*2d0**31, int64), 31), int(u(2)*2d0**31, int64))
         end do
      end do

      ex(1)%name = 'FILLING FACTOR';  ex(1)%note = '40×40 dish · 300 soups per row'
      ex(1)%nc = 19
      do i = 1, 19
         ex(1)%c(i) = cond_t(L=40, side=0, trials=300, fill=0.05d0*i)
      end do
      ex(2)%name = 'SIZE OF THE PETRI DISH';  ex(2)%note = 'fill 0.35 · 200 soups per row'
      ex(2)%nc = 9
      do i = 1, 9
         ex(2)%c(i) = cond_t(L=SIZES(i), side=0, trials=200, fill=0.35d0)
      end do
      ex(3)%name = 'LIFETIME DISTRIBUTION';  ex(3)%note = '40×40 dish · fill 0.35'
      ex(3)%nc = 1
      ex(3)%c(1) = cond_t(L=40, side=0, trials=NT3, fill=0.35d0)
      ex(4)%name = 'SEED COLONY SIZE';  ex(4)%note = 'fill 0.5 inside the colony · 40×40 dish'
      ex(4)%nc = 11
      do i = 1, 11
         ex(4)%c(i) = cond_t(L=40, side=i+1, trials=500, fill=0.5d0)
      end do
      do i = 1, 4
         ngames = ngames + sum(ex(i)%c(1:ex(i)%nc)%trials)
      end do
   end subroutine setup

   subroutine seed_rng(seed)
      ! spread one integer seed over the whole state of the intrinsic generator
      integer, intent(in) :: seed
      integer, allocatable :: put(:)
      integer :: n, i
      real(8) :: burn(16)
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = int(mod(int(seed, int64)*2654435761_int64 + 40503_int64*i, 2147483647_int64))
      end do
      call random_seed(put=put)
      call random_number(burn)
   end subroutine seed_rng

   ! ===================================================== the game itself

   subroutine new_game(g, cd, seed)
      ! pour a random soup (or a side x side seed colony) into an empty dish
      type(game_t), intent(inout) :: g
      type(cond_t), intent(in)    :: cd
      integer,      intent(in)    :: seed
      real(8), allocatable :: u(:,:)
      integer :: s, i0
      if (allocated(g%c)) then
         if (size(g%c, 1) /= cd%L + 2) deallocate(g%c, g%age)
      end if
      if (.not. allocated(g%c))   allocate(g%c(0:cd%L+1, 0:cd%L+1), g%age(cd%L, cd%L))
      if (.not. allocated(g%pop)) allocate(g%pop(0:GMAX))
      g%L = cd%L;  g%side = cd%side;  g%fill = cd%fill;  g%seed = seed
      g%gen = 0;  g%life = 0;  g%period = 0;  g%done = .false.

      s = cd%side
      if (s == 0) s = cd%L
      i0 = (cd%L - s)/2
      call seed_rng(seed)
      allocate(u(s, s))
      call random_number(u)
      g%c = 0
      g%c(i0+1:i0+s, i0+1:i0+s) = merge(1_int8, 0_int8, u < cd%fill)
      g%age = merge(1, 0, g%c(1:g%L, 1:g%L) == 1)
      g%pop(0) = count(g%c == 1)
      g%h(0) = iparity(zkey(1:g%L, 1:g%L), mask = g%c(1:g%L, 1:g%L) == 1)
   end subroutine new_game

   subroutine advance(g)
      ! one generation of the usual rules (birth on 3, survival on 2 or 3),
      ! then check whether the dish has stagnated
      type(game_t), intent(inout) :: g
      integer(int8)  :: nb(g%L, g%L)
      integer(int64) :: hash
      integer :: m, n, p
      m = g%L
      nb = g%c(0:m-1, 0:m-1) + g%c(1:m, 0:m-1) + g%c(2:m+1, 0:m-1) &
         + g%c(0:m-1, 1:m  )                   + g%c(2:m+1, 1:m  ) &
         + g%c(0:m-1, 2:m+1) + g%c(1:m, 2:m+1) + g%c(2:m+1, 2:m+1)
      g%c(1:m, 1:m) = merge(1_int8, 0_int8, nb == 3 .or. (nb == 2 .and. g%c(1:m, 1:m) == 1))
      g%age  = merge(g%age + 1, 0, g%c(1:m, 1:m) == 1)
      n      = count(g%c(1:m, 1:m) == 1)
      hash   = iparity(zkey(1:m, 1:m), mask = g%c(1:m, 1:m) == 1)
      g%gen  = g%gen + 1
      g%pop(g%gen) = n

      ! stagnant once a state recurs: compare with the last PMAX generations
      ! (64-bit Zobrist hash + population; a false match is ~1e-19 unlikely)
      do p = 1, min(g%gen, PMAX)
         if (g%h(mod(g%gen - p, PMAX)) == hash .and. g%pop(g%gen - p) == n) then
            g%done = .true.;  g%period = p;  g%life = g%gen - p
            exit
         end if
      end do
      g%h(mod(g%gen, PMAX)) = hash
      if (g%gen == GMAX .and. .not. g%done) then
         g%done = .true.;  g%period = 0;  g%life = GMAX
      end if
   end subroutine advance

   subroutine record(e, r, g)
      ! add a finished game to the tallies
      integer,      intent(in) :: e, r
      type(game_t), intent(in) :: g
      integer :: fin, k
      real(8) :: tl
      fin = g%pop(g%gen)
      tl  = g%life
      associate (t => ex(e)%t(r))
         t%n   = t%n + 1
         t%st  = t%st + tl
         t%st2 = t%st2 + tl*tl
         t%ash = t%ash + real(fin, 8)/real(g%L, 8)**2
         if (fin == 0)         t%extinct = t%extinct + 1
         if (g%period >= 2)    t%osc     = t%osc + 1
         if (g%period == 0)    t%stuck   = t%stuck + 1
         if (g%life >= TLONG)  t%long    = t%long + 1
         if (g%life > t%tmax) then
            t%tmax = g%life;  t%seedmax = g%seed
         end if
      end associate
      if (g%period == 0) then
         census(PMAX+1) = census(PMAX+1) + 1
      else if (fin == 0) then
         census(0) = census(0) + 1
      else
         census(g%period) = census(g%period) + 1
      end if
      if (e == 3) then                         ! insertion keeps life3 sorted
         k = n3
         do while (k > 0)
            if (life3(k) <= g%life) exit
            life3(k+1) = life3(k)
            k = k - 1
         end do
         life3(k+1) = g%life
         n3 = n3 + 1
      end if
   end subroutine record

   integer function next_seed()
      ! game k of a run gets seed (base*ngames + k) mod (2^31 - 1): different
      ! base seeds give disjoint sets of games, and any game can be replayed
      nplayed   = nplayed + 1
      next_seed = int(modulo(int(base, int64)*ngames + nplayed, 2147483647_int64))
   end function next_seed

   ! ===================================================== running it

   subroutine run_experiment(e)
      ! play every condition's games back to back as fast as possible; between
      ! chunks of work, show one of the batch's games evolving at watchable speed
      integer, intent(in) :: e
      type(game_t) :: bat, shw
      integer :: r, k, hold, f
      real(8) :: tnext, tstop
      if (.not. animate) write(output_unit, '(a)') '  experiment '//itoa(e)//' of 4: '// &
                                                   trim(ex(e)%name)
      r = 1;  k = 0;  hold = 0
      call new_game(bat, ex(e)%c(1), next_seed())
      if (animate) call new_game(shw, ex(e)%c(1), bat%seed)
      tnext = now()
      do while (r <= ex(e)%nc)
         tstop = now() + 0.75d0/FPS          ! the batch gets ~3/4 of every frame
         do
            call advance(bat)
            if (bat%done) then
               call record(e, r, bat)
               k = k + 1
               if (k == ex(e)%c(r)%trials) then
                  k = 0;  r = r + 1
                  if (r > ex(e)%nc) exit
               end if
               call new_game(bat, ex(e)%c(r), next_seed())
            end if
            if (animate .and. iand(bat%gen, 15) == 0) then
               if (now() > tstop) exit
            end if
         end do
         if (animate) call frame(e, r, shw, bat, hold, tnext)
      end do
      if (animate) then
         do f = 1, nint(3*FPS)               ! linger on the final tallies
            call frame(e, r, shw, bat, hold, tnext)
         end do
      end if
   end subroutine run_experiment

   subroutine frame(e, r, shw, bat, hold, tnext)
      ! one generation of the shown game, one redraw, then wait for the next tick
      integer,      intent(in)    :: e, r
      type(game_t), intent(inout) :: shw
      type(game_t), intent(in)    :: bat
      integer,      intent(inout) :: hold
      real(8),      intent(inout) :: tnext
      if (.not. shw%done) then
         call advance(shw)
      else
         hold = hold + 1                     ! pause on the outcome, then pick up
         if (hold > nint(1.5d0*FPS) .and. r <= ex(e)%nc) then   ! the batch's game
            call new_game(shw, ex(e)%c(r), bat%seed)
            hold = 0
         end if
      end if
      call draw(e, r, shw)
      tnext = max(tnext + 1/FPS, now() - 0.1d0)
      call sleep_until(tnext)
   end subroutine frame

   subroutine methuselah()
      ! replay the longest-lived soup of experiment 3, squeezed into ~12 s
      type(game_t) :: g
      real(8) :: tnext
      integer :: j, linger
      call new_game(g, ex(3)%c(1), ex(3)%t(1)%seedmax)
      mspeed = max(1, ceiling(ex(3)%t(1)%tmax/(12*FPS)))
      tnext = now()
      linger = 0
      do while (linger < nint(3*FPS))        ! run, then linger on the end
         if (g%done) then
            linger = linger + 1
         else
            do j = 1, mspeed
               if (.not. g%done) call advance(g)
            end do
         end if
         call draw(5, 0, g)
         tnext = max(tnext + 1/FPS, now() - 0.1d0)
         call sleep_until(tnext)
      end do
   end subroutine methuselah

   ! ===================================================== drawing

   subroutine draw(e, r, g)
      ! assemble the whole dashboard and write it in one go (no flicker)
      integer,      intent(in) :: e, r
      type(game_t), intent(in) :: g
      character(len=LW) :: lt(NROW), rt(NROW)
      integer :: i
      call dish_box(g, lt)
      call panel(e, r, g, rt)
      fbn = 0
      call put(ESC//'[?2026h'//ESC//'[H')     ! each line: erase it, then write it
      call put(ESC//'[K'//title_line(e)//NL)
      do i = 1, NROW
         call put(ESC//'[K'//trim(lt(i))//'  '//trim(rt(i))//NL)
      end do
      call put(ESC//'[K'//' '//pad(clip(status(g), CW + 2), CW + 3)//progress(e, g)//NL)
      call put(ESC//'[K'//'  '//pad(sparkline(g), CW + 2)//census_line()//NL)
      call put(ESC//'[J'//ESC//'[?2026l')
      call emit(fb(1:fbn))
   end subroutine draw

   subroutine dish_box(g, ln)
      ! the dish in half-block pixels: each character shows two cells stacked
      type(game_t),     intent(in)  :: g
      character(len=*), intent(out) :: ln(:)
      character(len=:), allocatable :: s, t
      integer :: f, off, row, col, top, bot, cf, cb
      f   = (g%L + CW - 1)/CW                ! big dishes are zoomed out by f
      off = (CW - (g%L + f - 1)/f)/2         ! small ones are centred
      t = 'PETRI DISH '//itoa(g%L)//'×'//itoa(g%L)
      if (f > 1) t = t//' at 1:'//itoa(f)
      ln(1) = box_top(t, CW)
      do row = 1, CW/2
         s = fg(C_RIM)//'│'
         cf = -9;  cb = -9                   ! colours in effect (-1: terminal's own)
         do col = 1, CW
            top = pixel(g, f, off, 2*row - 1, col)
            bot = pixel(g, f, off, 2*row,     col)
            if (top < 0 .and. bot < 0) then
               s = s//setc(cf, cb, cf, -1)//' '
            else if (top == bot) then        ! both halves alike: background only
               s = s//setc(cf, cb, cf, top)//' '
            else if (top < 0) then
               s = s//setc(cf, cb, bot, -1)//'▄'
            else
               s = s//setc(cf, cb, top, bot)//'▀'
            end if
         end do
         ln(row + 1) = s//RST//fg(C_RIM)//'│'//RST
      end do
      if (g%side > 0) then
         t = 'seed '//itoa(g%seed)//' · colony '//itoa(g%side)//'×'//itoa(g%side)
      else
         t = 'seed '//itoa(g%seed)//' · fill '//fmtf(g%fill, 4, 2)
      end if
      ln(NROW) = box_bottom(t, CW)
   end subroutine dish_box

   integer function pixel(g, f, off, pr, pc)
      ! colour of canvas pixel (pr,pc): -1 outside the dish, agar if empty,
      ! otherwise the colour of the youngest live cell it covers
      type(game_t), intent(in) :: g
      integer,      intent(in) :: f, off, pr, pc
      integer :: i0, i1, j0, j1, a
      if (min(pr, pc) <= off .or. max(pr, pc) > off + (g%L + f - 1)/f) then
         pixel = -1
         return
      end if
      i0 = (pr - off - 1)*f + 1;  i1 = min(i0 + f - 1, g%L)
      j0 = (pc - off - 1)*f + 1;  j1 = min(j0 + f - 1, g%L)
      a = minval(g%age(i0:i1, j0:j1), mask = g%age(i0:i1, j0:j1) > 0)
      if (a == huge(a)) then
         pixel = C_AGAR
      else
         pixel = C_AGE(1 + count(a > [1, 2, 4, 10, 40]))
      end if
   end function pixel

   function setc(cf, cb, nf, nbk) result(s)
      ! escape codes switching to foreground nf / background nbk, if not already set
      integer, intent(inout) :: cf, cb
      integer, intent(in)    :: nf, nbk
      character(len=:), allocatable :: s
      s = ''
      if (nf /= cf) then
         if (nf >= 0) s = fg(nf)
         cf = nf
      end if
      if (nbk /= cb) then
         if (nbk >= 0) then
            s = s//bgc(nbk)
         else
            s = s//ESC//'[49m'
         end if
         cb = nbk
      end if
   end function setc

   subroutine panel(e, r, g, ln)
      ! the box of tallies to the right of the dish
      integer,          intent(in)  :: e, r
      type(game_t),     intent(in)  :: g
      character(len=*), intent(out) :: ln(:)
      character(len=LW) :: inner(NROW - 2)
      character(len=:), allocatable :: title, foot
      integer :: i
      inner = ''
      select case (e)
      case (1, 2, 4)
         call table_panel(e, r, inner, foot)
      case (3)
         call hist_panel(inner, foot)
      case default
         call methuselah_panel(g, inner, foot)
      end select
      if (e <= 4) then
         title = itoa(e)//'/4 · '//trim(ex(e)%name)
      else
         title = 'THE METHUSELAH'
      end if
      ln(1) = box_top(title, PW)
      do i = 1, NROW - 2
         ln(i + 1) = fg(C_RIM)//'│'//RST//pad(clip(trim(inner(i)), PW), PW)//RST// &
                     fg(C_RIM)//'│'//RST
      end do
      ln(NROW) = box_bottom(foot, PW)
   end subroutine panel

   subroutine table_panel(e, r, ln, foot)
      ! one row per condition: mean lifetime (with a bar) and how the games ended
      integer,          intent(in)    :: e, r
      character(len=*), intent(inout) :: ln(:)
      character(len=:), allocatable, intent(out) :: foot
      character(len=:), allocatable :: s
      integer :: i, bw, best
      real(8) :: m, se, top, a

      bw = merge(12, 16, e == 4)
      top = 1;  best = 0
      do i = 1, ex(e)%nc
         if (ex(e)%t(i)%n == 0) cycle
         m = ex(e)%t(i)%st/ex(e)%t(i)%n
         if (m > top) then
            top = m;  best = i
         end if
      end do

      select case (e)
      case (1);      s = '  fill'
      case (2);      s = '     L'
      case default;  s = '  side cells'
      end select
      s = s//'  mean T ± s.e.  '//pad('lifetime', bw)//'  '
      if (e == 4) then
         s = s//lpad(itoa(TLONG)//'+', 4)//'  dead'
      else
         s = s//'dead     ash'
      end if
      ln(1) = fg(C_DIM)//s

      do i = 1, ex(e)%nc
         associate (t => ex(e)%t(i), c => ex(e)%c(i))
            select case (e)
            case (1);  s = fmtf(c%fill, 4, 2)
            case (2);  s = lpad(itoa(c%L), 4)
            case (4);  s = lpad(itoa(c%side), 4)//' '//fmtf(c%fill*c%side**2, 5, 1)
            end select
            if (i == r) then
               s = fg(C_HOT)//'▸ '//RST//BOLD//s//RST
            else
               s = '  '//s
            end if
            if (t%n == 0) then
               ln(i + 1) = s//fg(C_RIM)//'       ·'
               cycle
            end if
            m  = t%st/t%n
            se = sqrt(max(0d0, t%st2/t%n - m*m)/max(1, t%n - 1))
            s = s//'  '//fmtf(m, 6, 1)//fg(C_DIM)//' ±'//fmtf(se, 5, 1)//RST//'  '// &
                fg(merge(C_HOT, C_BAR, i == best))//pad(hbar(m/top, bw), bw)//RST//'  '
            if (e == 4) then
               s = s//ipct(t%long, t%n)//'  '//ipct(t%extinct, t%n)
            else
               s = s//ipct(t%extinct, t%n)//'  '//fmtf(100*t%ash/t%n, 5, 2)//'%'
            end if
            ln(i + 1) = s
         end associate
      end do

      foot = trim(ex(e)%note)
      if (e == 2) then                       ! power law through the finished rows
         a = fit_slope()
         if (a > 0) then
            ln(13) = '  power-law fit for L ≥ 16:   mean T ~ L^'//BOLD//fmtf(a, 4, 2)//RST
            ln(14) = fg(C_DIM)//'  doubling the dish multiplies T by '//fmtf(2**a, 4, 2)
         end if
      else if (e == 4) then                  ! smallest colony that usually lasts
         do i = 1, ex(4)%nc
            if (ex(4)%t(i)%n < ex(4)%c(i)%trials) exit
            if (nint(100d0*ex(4)%t(i)%long/ex(4)%t(i)%n) >= 50) then   ! as shown above
               ln(15) = '  colonies first last '//itoa(TLONG)//'+ generations half the time'
               ln(16) = '  at side '//BOLD//itoa(ex(4)%c(i)%side)//RST//' (about '// &
                        itoa(nint(ex(4)%c(i)%fill*ex(4)%c(i)%side**2))//' live cells)'
               exit
            end if
         end do
      end if
   end subroutine table_panel

   real(8) function fit_slope()
      ! least-squares slope of log(mean T) against log(L), finished rows, L >= 16
      real(8) :: x, y, sx, sy, sxx, sxy
      integer :: i, k
      sx = 0;  sy = 0;  sxx = 0;  sxy = 0;  k = 0
      do i = 1, ex(2)%nc
         if (ex(2)%c(i)%L < 16 .or. ex(2)%t(i)%n < ex(2)%c(i)%trials) cycle
         x = log(real(ex(2)%c(i)%L, 8))
         y = log(max(1d-9, ex(2)%t(i)%st/ex(2)%t(i)%n))
         sx = sx + x;  sy = sy + y;  sxx = sxx + x*x;  sxy = sxy + x*y;  k = k + 1
      end do
      fit_slope = -1
      if (k >= 2) fit_slope = (k*sxy - sx*sy)/(k*sxx - sx*sx)
   end function fit_slope

   subroutine hist_panel(ln, foot)
      ! histogram of the lifetimes of experiment 3, plus its key numbers
      integer, parameter :: NB = 46, HH = 12
      character(len=*), intent(inout) :: ln(:)
      character(len=:), allocatable, intent(out) :: foot
      character(len=LW) :: bars(HH)
      integer :: cnt(NB), i, w, b
      real(8) :: m, sd
      foot = trim(ex(3)%note)
      if (n3 == 0) return
      w = ceiling(real(life3(n3) + 1, 8)/NB)     ! bin width: NB bins reach the longest
      if (w > 10) w = 5*((w + 4)/5)              ! game, rounded up to a multiple of 5
      cnt = 0
      do i = 1, n3
         b = min(NB, life3(i)/w + 1)
         cnt(b) = cnt(b) + 1
      end do
      call vbars(cnt, HH, bars)
      ln(1) = fg(C_DIM)//'  generations a soup lasts before it stagnates'
      do i = 1, HH
         ln(i + 2) = '  '//fg(C_BAR)//trim(bars(i))
      end do
      ln(HH + 3) = '  '//fg(C_RIM)//repeat('─', NB)
      ln(HH + 4) = fg(C_DIM)//'  0'//lpad(itoa(w*NB/2), NB/2 - 1)//lpad(itoa(w*NB), NB/2)//' gen'
      m  = ex(3)%t(1)%st/n3
      sd = sqrt(max(0d0, ex(3)%t(1)%st2/n3 - m*m))
      ln(HH + 6) = '  median '//BOLD//itoa(life3((n3 + 1)/2))//RST//'   mean '//BOLD// &
                   itoa(nint(m))//RST//'   s.d. '//itoa(nint(sd))//'   soups '//itoa(n3)
      ln(HH + 7) = '  90% settle by gen '//itoa(life3(max(1, ceiling(0.9d0*n3))))
      ln(HH + 8) = '  longest '//BOLD//itoa(life3(n3))//RST//' gen'//fg(C_DIM)//', seed '// &
                   itoa(ex(3)%t(1)%seedmax)
   end subroutine hist_panel

   subroutine methuselah_panel(g, ln, foot)
      ! the longest-lived soup again, its population curve drawn as it goes
      integer, parameter :: NB = 46, HH = 11
      type(game_t),     intent(in)    :: g
      character(len=*), intent(inout) :: ln(:)
      character(len=:), allocatable, intent(out) :: foot
      character(len=LW) :: bars(HH)
      integer :: v(NB), j, a, b, span
      span = ex(3)%t(1)%tmax + 1
      ln(2) = '  the longest-lived of the '//itoa(n3)//' soups of experiment 3:'
      ln(3) = '  seed '//BOLD//itoa(g%seed)//RST//', active for '//BOLD//itoa(span - 1)//RST// &
              ' generations'
      ln(5) = fg(C_DIM)//'  population'
      v = 0
      do j = 1, NB
         a = (j - 1)*span/NB
         b = min(j*span/NB - 1, g%gen)
         if (b >= a) v(j) = maxval(g%pop(a:b))
      end do
      call vbars(v, HH, bars)
      do j = 1, HH
         ln(j + 5) = '  '//fg(C_BAR)//trim(bars(j))
      end do
      ln(HH + 6) = '  '//fg(C_RIM)//repeat('─', NB)
      ln(HH + 7) = fg(C_DIM)//'  0'//lpad(itoa(span - 1), NB - 1)//' gen'
      ln(HH + 9) = '  that is '//BOLD//fmtf(real(span - 1, 8)/max(1, life3((n3 + 1)/2)), 3, 1)// &
                   RST//' times the median soup ('//itoa(life3((n3 + 1)/2))//' gen)'
      foot = 'replay: '//itoa(mspeed)//' generations per frame'
   end subroutine methuselah_panel

   function title_line(e) result(s)
      integer, intent(in) :: e
      character(len=:), allocatable :: s, rhs
      integer :: secs
      secs = int(now() - tstart)
      if (e <= 4) then
         rhs = 'experiment '//itoa(e)//' of 4'
      else
         rhs = 'finale'
      end if
      rhs = rhs//' · '//itoa(secs/60)//':'//two(mod(secs, 60))//' '
      s = BOLD//fg(C_HOT)//' REINCARNATION'//RST//fg(C_DIM)// &
          "  ·  Conway's Game of Life, over and over"//RST
      s = pad(s, TW - vlen(rhs))//fg(C_DIM)//rhs//RST
   end function title_line

   function status(g) result(s)
      ! one line about the game in the dish
      type(game_t), intent(in) :: g
      character(len=:), allocatable :: s
      if (g%done .and. g%period == 0) then
         s = fg(C_HOT)//'no repeat within '//itoa(GMAX)//' generations'//RST
      else if (g%done) then
         s = fg(C_HOT)//'settled at gen '//itoa(g%life)//': '//fate(g)//RST
      else
         s = fg(C_DIM)//'gen '//RST//lpad(itoa(g%gen), 5)//fg(C_DIM)//'   population '//RST// &
             itoa(g%pop(g%gen))//fg(C_DIM)//' ('// &
             trim(adjustl(fmtf(100d0*g%pop(g%gen)/real(g%L, 8)**2, 4, 1)))//'%)'//RST
      end if
   end function status

   function fate(g) result(s)
      type(game_t), intent(in) :: g
      character(len=:), allocatable :: s
      if (g%period == 0) then
         s = 'still going at the cap'
      else if (g%pop(g%gen) == 0) then
         s = 'died out'
      else if (g%period == 1) then
         s = 'froze into still lifes'
      else
         s = 'oscillates, period '//itoa(g%period)
      end if
   end function fate

   function sparkline(g) result(s)
      ! population history of the game in the dish, squeezed into CW columns
      type(game_t), intent(in) :: g
      character(len=:), allocatable :: s
      character(len=LW) :: one(1)
      integer :: v(CW), j, ng
      ng = g%gen + 1
      v = 0
      do j = 1, CW
         if (ng <= CW) then
            if (j <= ng) v(j) = g%pop(j - 1)
         else
            v(j) = maxval(g%pop((j - 1)*ng/CW : j*ng/CW - 1))
         end if
      end do
      call vbars(v, 1, one)
      s = fg(C_BAR)//trim(one(1))//RST
   end function sparkline

   function progress(e, g) result(s)
      integer,      intent(in) :: e
      type(game_t), intent(in) :: g
      character(len=:), allocatable :: s, tag
      integer :: done, total
      if (e <= 4) then
         done  = sum(ex(e)%t(1:ex(e)%nc)%n)
         total = sum(ex(e)%c(1:ex(e)%nc)%trials)
         tag   = ' games'
      else
         done  = min(g%gen, ex(3)%t(1)%tmax)
         total = max(1, ex(3)%t(1)%tmax)
         tag   = ' gen'
      end if
      s = bgc(C_AGAR)//fg(C_BAR)//pad(hbar(real(done, 8)/total, 30), 30)//RST// &
          fg(C_DIM)//'  '//itoa(done)//' / '//itoa(total)//tag//RST
   end function progress

   function census_line() result(s)
      ! how all games so far have ended
      character(len=:), allocatable :: s
      integer :: tot
      tot = max(1, sum(census))
      s = fg(C_DIM)//'ended '//RST//'dead'//ipct(census(0), tot)//'  still'// &
          ipct(census(1), tot)//'  period 2'//ipct(census(2), tot)//'  other'// &
          ipct(sum(census(3:)), tot)
   end function census_line

   subroutine report()
      ! the final tallies, printed so they stay in the terminal
      character(len=LW) :: ln(NROW)
      character(len=:), allocatable :: s, blank
      type(game_t) :: none
      integer :: i, k, tot, last
      if (animate) call emit(ESC//'[2J'//ESC//'[H')
      tot = max(1, sum(census))
      write(output_unit, '(/,a,/)') BOLD//fg(C_HOT)//' REINCARNATION'//RST//fg(C_DIM)// &
           '  ·  results for base seed '//itoa(base)//'  ·  '//itoa(sum(census))// &
           ' games in '//itoa(nint(now() - tstart))//' s'//RST
      blank = fg(C_RIM)//'│'//RST//repeat(' ', PW)//RST//fg(C_RIM)//'│'//RST
      do k = 1, 4
         call panel(k, 0, none, ln)
         last = NROW - 1                       ! drop empty rows at the bottom of a box
         do while (last > 2 .and. ln(last) == blank)
            last = last - 1
         end do
         do i = 1, NROW
            if (i <= last + 1 .or. i == NROW) write(output_unit, '(a)') ' '//trim(ln(i))
         end do
      end do
      write(output_unit, '(/,a)') BOLD//'  how the '//itoa(sum(census))//' games ended'//RST
      do k = 0, PMAX + 1
         if (census(k) == 0) cycle
         select case (k)
         case (0);         s = 'died out'
         case (1);         s = 'still life'
         case (PMAX + 1);  s = 'unsettled'
         case default;     s = 'period '//itoa(k)
         end select
         write(output_unit, '(a)') '    '//pad(s, 12)//lpad(itoa(census(k)), 7)// &
              fmtf(100d0*census(k)/tot, 7, 2)//'%  '//fg(C_BAR)// &
              hbar(real(census(k), 8)/tot, 30)//RST
      end do
      write(output_unit, '(a)') ''
   end subroutine report

   ! ===================================================== small helpers

   subroutine vbars(v, h, ln)
      ! vertical bar chart of v, h text rows tall, in steps of 1/8 of a row
      integer,          intent(in)  :: v(:), h
      character(len=*), intent(out) :: ln(:)
      character(len=:), allocatable :: s
      integer :: i, k, lev, e8, vmax
      vmax = max(1, maxval(v))
      do k = 1, h
         lev = 8*(h - k)                     ! eighths below this row
         s = ''
         do i = 1, size(v)
            e8 = nint(8d0*h*v(i)/vmax)
            if (v(i) > 0) e8 = max(e8, 1)
            if (e8 >= lev + 8) then
               s = s//'█'
            else if (e8 > lev) then
               s = s//V8(e8 - lev)
            else
               s = s//' '
            end if
         end do
         ln(k) = s
      end do
   end subroutine vbars

   function hbar(x, w) result(s)
      ! horizontal bar for 0 <= x <= 1, w characters at full length
      real(8), intent(in) :: x
      integer, intent(in) :: w
      character(len=:), allocatable :: s
      integer :: n8
      n8 = nint(max(0d0, min(1d0, x))*w*8)
      s = repeat('█', n8/8)
      if (mod(n8, 8) > 0) s = s//H8(mod(n8, 8))
   end function hbar

   pure integer function vlen(s)
      ! width of s on screen: skip ESC[...m codes, count UTF-8 lead bytes only
      character(*), intent(in) :: s
      integer :: i
      logical :: in_esc
      vlen = 0;  in_esc = .false.
      do i = 1, len(s)
         if (in_esc) then
            in_esc = s(i:i) /= 'm'
         else if (s(i:i) == ESC) then
            in_esc = .true.
         else if (iand(iachar(s(i:i)), 192) /= 128) then
            vlen = vlen + 1
         end if
      end do
   end function vlen

   function pad(s, w) result(r)
      character(*), intent(in) :: s
      integer,      intent(in) :: w
      character(len=:), allocatable :: r
      r = s//repeat(' ', max(0, w - vlen(s)))
   end function pad

   function lpad(s, w) result(r)
      character(*), intent(in) :: s
      integer,      intent(in) :: w
      character(len=:), allocatable :: r
      r = repeat(' ', max(0, w - vlen(s)))//s
   end function lpad

   function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=16) :: buf
      write(buf, '(i0)') i
      s = trim(buf)
   end function itoa

   function two(i) result(s)
      integer, intent(in) :: i
      character(len=2) :: s
      write(s, '(i2.2)') i
   end function two

   function fmtf(x, w, d) result(s)
      ! x in a field of width w with d decimals (as an integer if it won't fit)
      real(8), intent(in) :: x
      integer, intent(in) :: w, d
      character(len=:), allocatable :: s
      character(len=32) :: buf, fmt
      if (abs(x) < 10d0**(w - d - 1) - 0.5d0*10d0**(-d)) then
         write(fmt, '(a,i0,a,i0,a)') '(f', w, '.', d, ')'
         write(buf, fmt) x
      else
         write(fmt, '(a,i0,a)') '(i', w, ')'
         write(buf, fmt) nint(x)
      end if
      s = buf(1:w)
   end function fmtf

   function ipct(k, n) result(s)
      integer, intent(in) :: k, n
      character(len=:), allocatable :: s
      s = lpad(itoa(nint(100d0*k/max(1, n))), 3)//'%'
   end function ipct

   function fg(c) result(s)                  ! 256-colour foreground
      integer, intent(in) :: c
      character(len=:), allocatable :: s
      s = ''
      if (color) s = ESC//'[38;5;'//itoa(c)//'m'
   end function fg

   function bgc(c) result(s)                 ! 256-colour background
      integer, intent(in) :: c
      character(len=:), allocatable :: s
      s = ''
      if (color) s = ESC//'[48;5;'//itoa(c)//'m'
   end function bgc

   function clip(s, w) result(r)
      ! s cut after w characters on screen, so a long label can never wrap a line
      character(*), intent(in) :: s
      integer,      intent(in) :: w
      character(len=:), allocatable :: r
      integer :: i, n
      logical :: in_esc
      n = 0;  in_esc = .false.
      do i = 1, len(s)
         if (in_esc) then
            in_esc = s(i:i) /= 'm'
         else if (s(i:i) == ESC) then
            in_esc = .true.
         else if (iand(iachar(s(i:i)), 192) /= 128) then
            n = n + 1
            if (n > w) then
               r = s(1:i-1)//RST
               return
            end if
         end if
      end do
      r = s
   end function clip

   function box_top(title, w) result(s)
      ! top edge of a box with inner width w, the title set into it
      character(*), intent(in) :: title
      integer,      intent(in) :: w
      character(len=:), allocatable :: s
      s = fg(C_RIM)//'╭─ '//BOLD//fg(C_TITLE)//title//RST//fg(C_RIM)//' '// &
          repeat('─', max(0, w - 3 - vlen(title)))//'╮'//RST
   end function box_top

   function box_bottom(label, w) result(s)
      ! bottom edge with a right-aligned label (if any)
      character(*), intent(in) :: label
      integer,      intent(in) :: w
      character(len=:), allocatable :: s
      if (len(label) == 0) then
         s = fg(C_RIM)//'╰'//repeat('─', w)//'╯'//RST
      else
         s = fg(C_RIM)//'╰'//repeat('─', max(0, w - 3 - vlen(label)))//' '//RST// &
             fg(C_DIM)//label//RST//fg(C_RIM)//' ─╯'//RST
      end if
   end function box_bottom

   subroutine put(s)
      character(*), intent(in) :: s
      fb(fbn + 1:fbn + len(s)) = s
      fbn = fbn + len(s)
   end subroutine put

   subroutine emit(s)
      character(*), intent(in) :: s
      write(output_unit, '(a)', advance='no') s
      flush(output_unit)
   end subroutine emit

   real(8) function now()
      integer(int64) :: cnt, rate
      call system_clock(cnt, rate)
      now = real(cnt, 8)/real(rate, 8)
   end function now

   subroutine sleep_until(t)
      real(8), intent(in) :: t
      real(8) :: dt
      integer(c_int) :: ierr
      dt = t - now()                         ! one clock reading: a negative argument
      if (dt > 0) ierr = usleep(int(min(dt, 1d0)*1d6, c_int))   ! would wrap around
   end subroutine sleep_until

end program reincarnation
