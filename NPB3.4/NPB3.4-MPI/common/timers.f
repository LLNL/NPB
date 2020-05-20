c---------------------------------------------------------------------
c---------------------------------------------------------------------
      
      module timers

      double precision start(64), elapsed(64)

      end module timers


c---------------------------------------------------------------------
c---------------------------------------------------------------------
      
      subroutine timer_clear(n)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      use timers
      implicit none

      integer n

      elapsed(n) = 0.0
      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine timer_start(n)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      use timers
      implicit none

      integer n

      include 'mpif.h'

      start(n) = MPI_Wtime()

      return
      end
      

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine timer_stop(n)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      use timers
      implicit none

      integer n

      include 'mpif.h'

      double precision t, now

      now = MPI_Wtime()
      t = now - start(n)
      elapsed(n) = elapsed(n) + t

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      double precision function timer_read(n)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      use timers
      implicit none

      integer n
      
      timer_read = elapsed(n)

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine check_timer_flag( timeron )

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none
      logical timeron

      integer nc, ios
      character(len=20) val

      timeron = .false.

c ... Check environment variable "NPB_TIMER_FLAG"
      call get_environment_variable('NPB_TIMER_FLAG', val, nc, ios)
      if (ios .eq. 0) then
         if (nc .le. 0) then
            timeron = .true.
         else if (val(1:1) .ge. '1' .and. val(1:1) .le. '9') then
            timeron = .true.
         else if (val .eq. 'on' .or. val .eq. 'ON' .or.
     &            val .eq. 'yes' .or. val .eq. 'YES' .or.
     &            val .eq. 'true' .or. val .eq. 'TRUE') then
            timeron = .true.
         endif

      else

c ... Check if the "timer.flag" file exists
         open (unit=2, file='timer.flag', status='old', iostat=ios)
         if (ios .eq. 0) then
            close(2)
            timeron = .true.
         endif

      endif

      return
      end

