!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.4         !
!                                                                         !
!                       O p e n M P     V E R S I O N                     !
!                                                                         !
!                                   L U                                   !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is an OpenMP version of the NPB LU code.              !
!    It is described in NAS Technical Report 99-011.                      !
!                                                                         !
!    Permission to use, copy, distribute and modify this software         !
!    for any purpose with or without fee is hereby granted.  We           !
!    request, however, that all derived work reference the NAS            !
!    Parallel Benchmarks 3.4. This software is provided "as is"           !
!    without express or implied warranty.                                 !
!                                                                         !
!    Information on NPB 3.4, including the technical report, the          !
!    original specifications, source code, results and information        !
!    on how to submit new results, is available at:                       !
!                                                                         !
!           http://www.nas.nasa.gov/Software/NPB/                         !
!                                                                         !
!    Send comments or suggestions to  npb@nas.nasa.gov                    !
!                                                                         !
!          NAS Parallel Benchmarks Group                                  !
!          NASA Ames Research Center                                      !
!          Mail Stop: T27A-1                                              !
!          Moffett Field, CA   94035-1000                                 !
!                                                                         !
!          E-mail:  npb@nas.nasa.gov                                      !
!          Fax:     (650) 604-3957                                        !
!                                                                         !
!-------------------------------------------------------------------------!

c---------------------------------------------------------------------
c
c Authors: S. Weeratunga
c          V. Venkatakrishnan
c          E. Barszcz
c          M. Yarrow
c          H. Jin
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------
      program applu
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c   driver for the performance evaluation of the solver for
c   five coupled parabolic/elliptic partial differential equations.
c
c---------------------------------------------------------------------

      use lu_data
      implicit none

      character class
      logical verified
      double precision mflops

      double precision t, tmax, timer_read, trecs(t_last)
      external timer_read
      integer i
      character t_names(t_last)*8

c---------------------------------------------------------------------
c     Setup info for timers
c---------------------------------------------------------------------

      call check_timer_flag( timeron )
      if (timeron) then
         t_names(t_total) = 'total'
         t_names(t_rhsx) = 'rhsx'
         t_names(t_rhsy) = 'rhsy'
         t_names(t_rhsz) = 'rhsz'
         t_names(t_rhs) = 'rhs'
         t_names(t_jacld) = 'jacld'
         t_names(t_blts) = 'blts'
         t_names(t_jacu) = 'jacu'
         t_names(t_buts) = 'buts'
         t_names(t_add) = 'add'
         t_names(t_l2norm) = 'l2norm'
      endif

c---------------------------------------------------------------------
c   read input data
c---------------------------------------------------------------------
      call read_input()

c---------------------------------------------------------------------
c   set up domain sizes
c---------------------------------------------------------------------
      call domain()

      call alloc_space

c---------------------------------------------------------------------
c   set up coefficients
c---------------------------------------------------------------------
      call setcoeff()

c---------------------------------------------------------------------
c   set the boundary values for dependent variables
c---------------------------------------------------------------------
      call setbv()

c---------------------------------------------------------------------
c   set the initial values for dependent variables
c---------------------------------------------------------------------
      call setiv()

c---------------------------------------------------------------------
c   compute the forcing term based on prescribed exact solution
c---------------------------------------------------------------------
      call erhs()

c---------------------------------------------------------------------
c   perform one SSOR iteration to touch all data pages
c---------------------------------------------------------------------
      call ssor(1)

c---------------------------------------------------------------------
c   reset the boundary and initial values
c---------------------------------------------------------------------
      call setbv()
      call setiv()

c---------------------------------------------------------------------
c   perform the SSOR iterations
c---------------------------------------------------------------------
      call ssor(itmax)

c---------------------------------------------------------------------
c   compute the solution error
c---------------------------------------------------------------------
      call error()

c---------------------------------------------------------------------
c   compute the surface integral
c---------------------------------------------------------------------
      call pintgr()

c---------------------------------------------------------------------
c   verification test
c---------------------------------------------------------------------
      call verify ( rsdnm, errnm, frc, class, verified )
      mflops = 1.0d-6*dble(itmax)*(1984.77*dble( nx0 )
     >     *dble( ny0 )
     >     *dble( nz0 )
     >     -10923.3*(dble( nx0+ny0+nz0 )/3.)**2 
     >     +27770.9* dble( nx0+ny0+nz0 )/3.
     >     -144010.)
     >     / maxtime

      call print_results('LU', class, nx0,
     >  ny0, nz0, itmax,
     >  maxtime, mflops, '          floating point', verified, 
     >  npbversion, compiletime, cs1, cs2, cs3, cs4, cs5, cs6, 
     >  '(none)')

c---------------------------------------------------------------------
c      More timers
c---------------------------------------------------------------------
      if (.not.timeron) goto 999

      do i=1, t_last
         trecs(i) = timer_read(i)
      end do
      tmax = maxtime
      if ( tmax .eq. 0. ) tmax = 1.0

      write(*,800)
 800  format('  SECTION     Time (secs)')
      do i=1, t_last
         if (i.ne.t_jacld .and. i.ne.t_jacu) then
            write(*,810) t_names(i), trecs(i), trecs(i)*100./tmax
         endif
         if (i.eq.t_rhs) then
            t = trecs(t_rhsx) + trecs(t_rhsy) + trecs(t_rhsz)
            write(*,820) 'sub-rhs', t, t*100./tmax
            t = trecs(i) - t
            write(*,820) 'rest-rhs', t, t*100./tmax
         endif
 810     format(2x,a8,':',f9.3,'  (',f6.2,'%)')
 820     format(5x,'--> ',a8,':',f9.3,'  (',f6.2,'%)')
      end do

 999  continue
      end


