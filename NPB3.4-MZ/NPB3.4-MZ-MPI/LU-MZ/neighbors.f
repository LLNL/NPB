c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine neighbors ()

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      use lu_data
      use mpinpb

      implicit none

c---------------------------------------------------------------------
c     figure out the neighbors and their wrap numbers for each processor
c---------------------------------------------------------------------

      south = -1
      east  = -1
      north = -1
      west  = -1

      if (row.gt.0) then
          north = myid - 1
      else
          north = -1
      end if

      if (row.lt.xdim-1) then
          south = myid + 1
      else
          south = -1
      end if

      if (col.gt.0) then
          west = myid - xdim
      else
          west = -1
      end if

      if (col.lt.zdim-1) then
          east = myid + xdim
      else 
          east = -1
      end if

      return
      end
