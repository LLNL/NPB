c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  work_lhs module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module work_lhs

      use sp_data, only : problem_size

c-----------------------------------------------------------------------
c   Working array for LHS
c-----------------------------------------------------------------------

      include 'blk_par.h'

      integer, parameter :: IMAXP=problem_size/2*2
      double precision
     >      lhs (blkdim,5,0:IMAXP),
     >      lhsp(blkdim,5,0:IMAXP),
     >      lhsm(blkdim,5,0:IMAXP),
     >      rhsx(blkdim,5,0:IMAXP),
     >      cv  (blkdim,0:problem_size-1),
     >      rhov(blkdim,0:problem_size-1)
!$omp threadprivate(lhs, lhsp, lhsm, rhsx, cv, rhov)

      end module work_lhs


c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine lhsinit(ni)

       use work_lhs
       implicit none

       integer ni

       integer i, j

c---------------------------------------------------------------------
c     zap the whole left hand side for starters
c     set all diagonal values to 1. This is overkill, but convenient
c---------------------------------------------------------------------
       do i = 0, ni, ni
          do j = 1, bsize
             lhs (j,1,i) = 0.0d0
             lhs (j,2,i) = 0.0d0
             lhs (j,3,i) = 1.0d0
             lhs (j,4,i) = 0.0d0
             lhs (j,5,i) = 0.0d0
             lhsp(j,1,i) = 0.0d0
             lhsp(j,2,i) = 0.0d0
             lhsp(j,3,i) = 1.0d0
             lhsp(j,4,i) = 0.0d0
             lhsp(j,5,i) = 0.0d0
             lhsm(j,1,i) = 0.0d0
             lhsm(j,2,i) = 0.0d0
             lhsm(j,3,i) = 1.0d0
             lhsm(j,4,i) = 0.0d0
             lhsm(j,5,i) = 0.0d0
          end do
       end do
 
       return
       end

