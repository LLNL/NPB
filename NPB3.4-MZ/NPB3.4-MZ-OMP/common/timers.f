c---------------------------------------------------------------------
c---------------------------------------------------------------------
      
      module timers

      double precision start(64), elapsed(64)
c$omp threadprivate(start, elapsed)

      double precision, external :: elapsed_time

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

      start(n) = elapsed_time()

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

      double precision t, now

      now = elapsed_time()
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

      double precision function elapsed_time()

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none
c$    external         omp_get_wtime
c$    double precision omp_get_wtime

      double precision t
      logical          mp

c ... Use the OpenMP timer if we can (via C$ conditional compilation)
      mp = .false.
c$    mp = .true.
c$    t = omp_get_wtime()

      if (.not.mp) then
c This function must measure wall clock time, not CPU time. 
c Since there is no portable timer in Fortran (77)
c we call a routine compiled in C (though the C source may have
c to be tweaked). 
         call wtime(t)
c The following is not ok for "official" results because it reports
c CPU time not wall clock time. It may be useful for developing/testing
c on timeshared Crays, though. 
c        call second(t)
      endif

      elapsed_time = t

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine check_timer_flag( timeron )

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none
      integer timeron

      integer ios
      character(len=20) val

      timeron = 0

c ... Check environment variable "NPB_TIMER_FLAG"
      call get_menv( 'NPB_TIMER_FLAG', val, ios )
      if (ios .ge. 0) then
         if (ios .eq. 0) then
            timeron = 1
         else if (val(1:1) .ge. '1' .and. val(1:1) .le. '9') then
            read (val,*,iostat=ios) timeron
            if (ios .ne. 0) timeron = 1
         else if (val .eq. 'on' .or. val .eq. 'ON' .or.
     &            val .eq. 'yes' .or. val .eq. 'YES' .or.
     &            val .eq. 'true' .or. val .eq. 'TRUE') then
            timeron = 1
         endif

      else

c ... Check if the "timer.flag" file exists
         open (unit=2, file='timer.flag', status='old', iostat=ios)
         if (ios .eq. 0) then
            read (2,*,iostat=ios) timeron
            if (ios .ne. 0) timeron = 1
            close(2)
         endif

      endif

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine get_menv( env, val, nc )
c
c  Get value for an environment variable (using intrisic routine)
c
      implicit none
c
      character(len=*) env, val
      integer nc
c
      integer ios
c
      call get_environment_variable(env, val, nc, ios)
      if (ios .ne. 0) nc = -1        ! env not defined
c
      return
      end

