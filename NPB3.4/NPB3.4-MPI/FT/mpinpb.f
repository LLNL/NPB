c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  mpinpb module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module mpinpb

      include 'mpif.h'

c--------------------------------------------------------------------
c 'np' number of processors, 'np_min' min number of processors
c--------------------------------------------------------------------
      integer np_min, np

c we need a bunch of logic to keep track of how
c arrays are laid out. 
c coords of this processor
      integer me, me1, me2

c need a communicator for row/col in processor grid
      integer commslice1, commslice2

c mpi data types
      integer dc_type

      end module mpinpb

