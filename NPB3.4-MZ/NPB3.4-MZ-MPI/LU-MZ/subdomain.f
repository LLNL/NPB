
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine subdomain(nx1, ny1, nz1, nx)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      use lu_data
      use mpinpb

      implicit none

      integer nx1, ny1, nz1, nx

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer mm, errorcode


c---------------------------------------------------------------------
c
c   set up the sub-domain sizes
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   original dimensions
c---------------------------------------------------------------------
      nx0 = nx1
      ny0 = ny1
      nz0 = nz1

c---------------------------------------------------------------------
c   x dimension
c---------------------------------------------------------------------
      mm   = mod(nx0,xdim)
      if (row.lt.mm) then
        nx = nx0/xdim + 1
        ipt = row*nx
      else
        nx = nx0/xdim
        ipt = row*nx + mm
      end if

c---------------------------------------------------------------------
c   check the sub-domain size
c---------------------------------------------------------------------
      if ( nx .lt. 3 ) then
         if (myid .eq. root) write (*,2001) nx, ny0, nz0
 2001    format (5x,'SUBDOMAIN SIZE IS TOO SMALL - ',
     >        /5x,'ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS',
     >        /5x,'SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL',
     >        /5x,'TO 3 THEY ARE CURRENTLY', 3I5)
         call error_cond( 0, ' ' )
      end if


c---------------------------------------------------------------------
c   set up the start and end in i extent for all processors
c---------------------------------------------------------------------
      ist = 1
      iend = nx
      if (north.eq.-1) ist = 2
      if (south.eq.-1) iend = nx - 1

      return
      end


