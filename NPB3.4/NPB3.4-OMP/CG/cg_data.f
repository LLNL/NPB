c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  cg_data module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module cg_data

      include 'npbparams.h'

c---------------------------------------------------------------------
c  Class specific parameters are defined in the npbparams.h
c  include file, which is written by the sys/setparams.c program.
c---------------------------------------------------------------------


c ... dimension parameters
      integer    nz, naz
      parameter( nz = na*(nonzer+1)*(nonzer+1) )
      parameter( naz = na*(nonzer+1) )

c ... main_int_mem
      integer, allocatable ::  colidx(:), rowstr(:),
     >                         iv(:),  arow(:), acol(:)

c ... main_flt_mem
      double precision, allocatable ::
     >                         v(:), aelt(:), a(:),
     >                         x(:),
     >                         z(:),
     >                         p(:),
     >                         q(:),
     >                         r(:)

c ... partition size
      integer                  naa, nzz, 
     >                         firstrow, 
     >                         lastrow, 
     >                         firstcol, 
     >                         lastcol

      double precision         amult, tran
!$omp threadprivate (amult, tran)

      external         timer_read
      double precision timer_read

      integer T_init, T_bench, T_conj_grad, T_last
      parameter (T_init=1, T_bench=2, T_conj_grad=3, T_last=3)

      logical timeron

      end module cg_data


c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  tinfo module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module tinfo

      integer        max_threads
      parameter      (max_threads=1024)

      integer        last_n(0:max_threads)

      integer        myid, num_threads, ilow, ihigh
!$omp threadprivate (myid, num_threads, ilow, ihigh)

      end module tinfo


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine alloc_space

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c allocate space dynamically for data arrays
c---------------------------------------------------------------------

      use cg_data
      implicit none

      integer ios


      allocate ( 
     >          colidx(nz), rowstr(na+1),
     >          iv(nz+na),  arow(na), acol(naz),
     >          v(nz), aelt(naz), a(nz),
     >          x(na+2),
     >          z(na+2),
     >          p(na+2),
     >          q(na+2),
     >          r(na+2),
     >          stat = ios)

      if (ios .ne. 0) then
         write(*,*) 'Error encountered in allocating space'
         stop
      endif

      return
      end

