c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  lu_data module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module lu_data

c---------------------------------------------------------------------
c   npbparams.h defines parameters that depend on the class and 
c   number of nodes
c---------------------------------------------------------------------

      include 'npbparams.h'

      integer(kind2), parameter :: 
     &            max_xysize=(gx_size+x_zones)*gy_size,
     &            proc_max_size=max_xysize*gz_size,
     &            proc_max_size5=proc_max_size*5,
     &            proc_max_bcsize=max_xysize*20

c---------------------------------------------------------------------
c   parameters which can be overridden in runtime config file
c   isiz1,isiz2,isiz3 give the maximum size
c   ipr = 1 to print out verbose information
c   omega = 2.0 is correct for all classes
c   tolrsd is tolerance levels for steady state residuals
c---------------------------------------------------------------------
      integer ipr_default
      parameter (ipr_default = 1)
      double precision omega_default
      parameter (omega_default = 1.2d0)
      double precision tolrsd1_def, tolrsd2_def, tolrsd3_def, 
     >                 tolrsd4_def, tolrsd5_def
      parameter (tolrsd1_def=1.0e-08, 
     >          tolrsd2_def=1.0e-08, tolrsd3_def=1.0e-08, 
     >          tolrsd4_def=1.0e-08, tolrsd5_def=1.0e-08)

      double precision c1, c2, c3, c4, c5
      parameter( c1 = 1.40d+00, c2 = 0.40d+00,
     >           c3 = 1.00d-01, c4 = 1.00d+00,
     >           c5 = 1.40d+00 )

c---------------------------------------------------------------------
c   grid
c---------------------------------------------------------------------
      double precision  dxi, deta, dzeta
      double precision  tx1, tx2, tx3
      double precision  ty1, ty2, ty3
      double precision  tz1, tz2, tz3

c---------------------------------------------------------------------
c   dissipation
c---------------------------------------------------------------------
      double precision dx1, dx2, dx3, dx4, dx5
      double precision dy1, dy2, dy3, dy4, dy5
      double precision dz1, dz2, dz3, dz4, dz5
      double precision dssp

      integer   max_zones
      parameter (max_zones=x_zones*y_zones)
      integer   x_start(x_zones), x_end(x_zones), x_size(x_zones),
     >          y_start(y_zones), y_end(y_zones), y_size(y_zones),
     >          iz_west (max_zones), iz_east (max_zones),
     >          iz_south(max_zones), iz_north(max_zones)

      integer(kind2) :: start1(max_zones), start5(max_zones),
     $          qstart_west (max_zones), qstart_east (max_zones),
     $          qstart_south(max_zones), qstart_north(max_zones)

c---------------------------------------------------------------------
c   output control parameters
c---------------------------------------------------------------------
      integer ipr, inorm, npb_verbose

c---------------------------------------------------------------------
c   newton-raphson iteration control parameters
c---------------------------------------------------------------------
      integer itmax, invert
      double precision  dt, omega, tolrsd(5), ttotal

c---------------------------------------------------------------------
c   coefficients of the exact solution
c---------------------------------------------------------------------
      double precision ce(5,13)

c---------------------------------------------------------------------
c   1-d working arrays
c---------------------------------------------------------------------
      double precision  flux(5,problem_size), 
     >                  utmp(6,problem_size),
     >                  rtmp(5,problem_size)
!$omp threadprivate( flux, utmp, rtmp )

c---------------------------------------------------------------------
c   timers
c---------------------------------------------------------------------
      integer t_rhsx,t_rhsy,t_rhsz,t_rhs,t_jacld,t_blts,t_jacu,
     >        t_buts,t_add,t_l2norm,t_rdis1,t_rdis2,t_last,t_total
      parameter (t_total = 1)
      parameter (t_rhsx = 2)
      parameter (t_rhsy = 3)
      parameter (t_rhsz = 4)
      parameter (t_rhs = 5)
      parameter (t_jacld = 6)
      parameter (t_blts = 7)
      parameter (t_jacu = 8)
      parameter (t_buts = 9)
      parameter (t_add = 10)
      parameter (t_l2norm = 11)
      parameter (t_rdis1 = 12)
      parameter (t_rdis2 = 13)
      parameter (t_last = 13)
      logical timeron
      double precision maxtime

      end module lu_data


c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  lu_fields module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module lu_fields

      use lu_data, only : max_zones, problem_size

      integer   nx(max_zones), nxmax(max_zones), ny(max_zones), 
     $          nz(max_zones),
     &          proc_zone_id(max_zones), proc_num_zones

      integer   isync(0:problem_size)

c---------------------------------------------------------------------
c   Define all field arrays as one-dimenional arrays, to be reshaped
c---------------------------------------------------------------------

      double precision , allocatable ::
     >                 u     (:),
     >                 rsd   (:),
     >                 frct  (:),
     >                 qs    (:),
     >                 rho_i (:),
     >                 qbc   (:)

c---------------------------------------------------------------------
c   2D auxiliary arrays are dimensioned to accommodate the largest
c   zone cross section
c---------------------------------------------------------------------

      double precision, allocatable ::
     $                 a   (:), 
     $                 b   (:), 
     $                 c   (:), 
     $                 d   (:),
     $                 au  (:), 
     $                 bu  (:), 
     $                 cu  (:), 
     $                 du  (:),
     $                 tv  (:),
     $                 phi1(:),
     $                 phi2(:)

      end module lu_fields


c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  Allocate space for field arrays
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine alloc_field_space

      use lu_data, only :
     &            proc_max_size,
     &            proc_max_size5,
     &            proc_max_bcsize
      use lu_fields

      implicit none

      integer ios

      allocate (
     >          u     (proc_max_size5),
     >          rsd   (proc_max_size5),
     >          frct  (proc_max_size5),
     >          qs    (proc_max_size ),
     >          rho_i (proc_max_size ),
     >          qbc   (proc_max_bcsize),
     >          stat=ios )

      if (ios .eq. 0) allocate (
     $          a (25*problem_size*problem_size), 
     $          b (25*problem_size*problem_size), 
     $          c (25*problem_size*problem_size), 
     $          d (25*problem_size*problem_size),
     $          au(25*problem_size*problem_size), 
     $          bu(25*problem_size*problem_size), 
     $          cu(25*problem_size*problem_size), 
     $          du(25*problem_size*problem_size),
     $          tv( 5*problem_size*problem_size),
     $          phi1 (problem_size*problem_size),
     $          phi2 (problem_size*problem_size),
     >          stat=ios )

      if (ios .ne. 0) then
         write(*,*) 'Error encountered in allocating space'
         stop
      endif

      return
      end

