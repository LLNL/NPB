!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.4         !
!                                                                         !
!          O p e n M P    M U L T I - Z O N E    V E R S I O N            !
!                                                                         !
!                            B T - M Z - O M P                            !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is an OpenMP version of the NPB BT code.              !
!    Refer to NAS Technical Reports 95-020 and 99-011 for details.        !
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
c Authors: R. Van der Wijngaart
c          T. Harris
c          M. Yarrow
c          H. Jin
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       program BT_MZ
c---------------------------------------------------------------------

       use bt_data
       use bt_fields
       use ompnpb

       implicit none

       integer          num_zones

       integer          i, niter, step, fstatus, zone, iz, itimer,
     >                  tot_threads, nthreads
       double precision navg, mflops, nsur, n3

       external         timer_read
       double precision tmax, timer_read, t, trecs(t_last),
     >                  tsum(t_last), tming(t_last), tmaxg(t_last)
       logical          verified
       character        t_names(t_last)*8


c---------------------------------------------------------------------
c      Reads input file (if it exists) else takes
c      defaults from parameters
c---------------------------------------------------------------------

       write(*, 1000)

       call check_timer_flag( itimer )

       open (unit=2,file='inputbt-mz.data',status='old', 
     >       iostat=fstatus)

       if (fstatus .eq. 0) then
         write(*,*) 'Reading from input file inputbt-mz.data'
         read (2,*) niter
         read (2,*) dt
         read (2,*,err=20,end=20) itimer
   20    close(2)

         if (niter .eq. 0)  niter = niter_default
         if (dt .eq. 0.d0)  dt    = dt_default

       else
         niter = niter_default
         dt    = dt_default
       endif

       write(*, 1001) x_zones, y_zones
       write(*, 1002) gx_size, gy_size, gz_size
       write(*, 1003) niter, dt
 1000  format(//, ' NAS Parallel Benchmarks (NPB3.4-MZ OpenMP)',
     >            ' - BT-MZ Benchmark', /)
 1001  format(' Number of zones: ', i3, ' x ', i3)
 1002  format(' Total mesh size: ', i5, ' x ', i5, ' x ', i3)
 1003  format(' Iterations: ', i3, '    dt: ', F10.6/)

       timeron = (itimer .gt. 0)
       if (timeron) then
         t_names(t_total) = 'total'
         t_names(t_rhsx) = 'rhsx'
         t_names(t_rhsy) = 'rhsy'
         t_names(t_rhsz) = 'rhsz'
         t_names(t_rhs) = 'rhs'
         t_names(t_xsolve) = 'xsolve'
         t_names(t_ysolve) = 'ysolve'
         t_names(t_zsolve) = 'zsolve'
         t_names(t_rdis1) = 'qbc_copy'
         t_names(t_rdis2) = 'qbc_comm'
         t_names(t_add) = 'add'
       endif

       call env_setup(tot_threads)

       call zone_setup(nx, nxmax, ny, nz)

       num_zones = max_zones
       call setup_omp(num_zones, nx, ny, nz, tot_threads)
       call zone_starts(num_zones, nx, nxmax, ny, nz)

c---------------------------------------------------------------------
c      allocate space for field arrays
c---------------------------------------------------------------------
       call alloc_field_space

       call set_constants


       if (timeron) then
          do i = 1, t_last
             tsum(i) = 0.d0
             tming(i) = huge(0.d0)
             tmaxg(i) = 0.d0
          end do
       endif

c---------------------------------------------------------------------
c      start of the outer parallel region
c---------------------------------------------------------------------
!$omp parallel private(iz,i,zone,step,t,trecs,nthreads,
!$omp&  proc_num_zones,proc_zone_id)
!$omp&  if(nested.ne.2)

       call init_omp(num_zones, proc_zone_id, proc_num_zones)

       nthreads = proc_num_threads(myid+1)
c$     call omp_set_num_threads(nthreads)

       do iz = 1, proc_num_zones
         zone = proc_zone_id(iz)

         call initialize(u(start5(zone)),
     $                   nx(zone), nxmax(zone), ny(zone), nz(zone))
         call exact_rhs(forcing(start5(zone)),
     $                  nx(zone), nxmax(zone), ny(zone), nz(zone))

       end do

       do i = 1, t_last
          call timer_clear(i)
       end do

c---------------------------------------------------------------------
c      do one time step to touch all code, and reinitialize
c---------------------------------------------------------------------

       call exch_qbc(u, qbc, nx, nxmax, ny, nz, 
     &               proc_zone_id, proc_num_zones)

       do iz = 1, proc_num_zones
         zone = proc_zone_id(iz)
         call adi(rho_i(start1(zone)), us(start1(zone)), 
     $            vs(start1(zone)), ws(start1(zone)), 
     $            qs(start1(zone)), square(start1(zone)), 
     $            rhs(start5(zone)), forcing(start5(zone)), 
     $            u(start5(zone)), 
     $            nx(zone), nxmax(zone), ny(zone), nz(zone))
       end do

       do iz = 1, proc_num_zones
         zone = proc_zone_id(iz)
         call initialize(u(start5(zone)), 
     $                   nx(zone), nxmax(zone), ny(zone), nz(zone)) 
       end do

       do i = 1, t_last
          call timer_clear(i)
       end do
!$omp barrier
       call timer_start(1)

c---------------------------------------------------------------------
c      start the benchmark time step loop
c---------------------------------------------------------------------

       do  step = 1, niter

!$omp master
         if (mod(step, 20) .eq. 0 .or. step .eq. 1) then
            write(*, 200) step
 200        format(' Time step ', i4)
         endif
!$omp end master

         call exch_qbc(u, qbc, nx, nxmax, ny, nz, 
     &                 proc_zone_id, proc_num_zones)

         do iz = 1, proc_num_zones
           zone = proc_zone_id(iz)
           call adi(rho_i(start1(zone)), us(start1(zone)), 
     $              vs(start1(zone)), ws(start1(zone)), 
     $              qs(start1(zone)), square(start1(zone)), 
     $              rhs(start5(zone)), forcing(start5(zone)), 
     $              u(start5(zone)), 
     $              nx(zone), nxmax(zone), ny(zone), nz(zone))
         end do

       end do

!$omp barrier
       call timer_stop(1)
       t = timer_read(1)

c---------------------------------------------------------------------
c      perform verification and print results
c---------------------------------------------------------------------

       call verify(niter, verified, num_zones, rho_i, us, vs, ws, 
     $             qs, square, rhs, forcing, u, nx, nxmax, ny, nz, 
     $             proc_zone_id, proc_num_zones)

!$omp master
       tmax = t
       mflops = 0.0d0
       if( tmax .ne. 0. ) then
         do zone = 1, num_zones
           n3 = dble(nx(zone))*ny(zone)*nz(zone)
           navg = (nx(zone) + ny(zone) + nz(zone))/3.0
           nsur = (nx(zone)*ny(zone) + nx(zone)*nz(zone) +
     >             ny(zone)*nz(zone))/3.0
           mflops = mflops + 1.0d-6*float(niter) *
     >      (3478.8d0 * n3 - 17655.7d0 * nsur + 28023.7d0 * navg)
     >      / tmax
         end do
       endif

       call print_results('BT-MZ', class, gx_size, gy_size, gz_size, 
     >                    niter, tmax, mflops, num_othreads, 
     >                    tot_threads, '          floating point', 
     >                    verified, npbversion,compiletime, cs1, cs2, 
     >                    cs3, cs4, cs5, cs6, '(none)')
!$omp end master
!$omp barrier

c---------------------------------------------------------------------
c      More timers
c---------------------------------------------------------------------
       if (.not.timeron) goto 999

       do i=1, t_last
          trecs(i) = timer_read(i)
       end do
       if (tmax .eq. 0.0) tmax = 1.0

       do i=1, t_last
!$omp atomic
         tsum(i) = tsum(i) + trecs(i)
!$omp atomic
         tming(i) = min(tming(i), trecs(i))
!$omp atomic
         tmaxg(i) = max(tmaxg(i), trecs(i))
       end do
!$omp barrier

!$omp master
!$     write(*, 700) num_othreads
!$     do i = 1, t_last
!$        tsum(i) = tsum(i) / num_othreads
!$        write(*, 710) i, t_names(i), tming(i), tmaxg(i), tsum(i)
!$     end do
!$omp end master
 700   format(' #othrs =', i6, 11x, 'minimum', 5x, 'maximum', 
     >        5x, 'average')
 710   format(' timer ', i2, '(', A8, ') :', 3(2x,f10.4))

!$     if (itimer .lt. 2) goto 999

!$omp barrier
!$omp critical (ptime)
!$     write(*,800) myid, nthreads
 800   format(/' myid =',i5,'   num_ithreads =',i4)
       write(*,805)
 805   format('  SECTION   Time (secs)')
       do i=1, t_last
          write(*,810) t_names(i), trecs(i), trecs(i)*100./tmax
          if (i.eq.t_rhs) then
             t = trecs(t_rhsx) + trecs(t_rhsy) + trecs(t_rhsz)
             write(*,820) 'sub-rhs', t, t*100./tmax
             t = trecs(t_rhs) - t
             write(*,820) 'rest-rhs', t, t*100./tmax
          elseif (i.eq.t_rdis2) then
             t = trecs(t_rdis1) + trecs(t_rdis2)
             write(*,820) 'exch_qbc', t, t*100./tmax
          endif
 810      format(2x,a8,':',f9.3,'  (',f6.2,'%)')
 820      format('    --> total ',a8,':',f9.3,'  (',f6.2,'%)')
       end do
!$omp end critical (ptime)

 999   continue

!$omp end parallel

       end

