c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  mpinpb module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module mpinpb

      include 'mpif.h'

      integer  node, no_nodes, total_nodes, root, comm_setup, 
     >         comm_solve, comm_rhs, dp_type
      logical  active

      end module mpinpb

