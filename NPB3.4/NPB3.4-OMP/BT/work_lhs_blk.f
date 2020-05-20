c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  work_lhs module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module work_lhs

      use bt_data, only : problem_size

      include 'blk_par.h'

      double precision fjac(blkdim, 5, 5, 0:2),
     >                 njac(blkdim, 5, 5, 0:2),
     >                 lhsa(blkdim, 5, 5, 0:2),
     >                 lhsb(blkdim, 5, 5, 0:2),
     >                 lhsc(blkdim, 5, 5, 0:problem_size-1),
     >                 rhsx(blkdim, 5, 0:problem_size-1)
!$omp threadprivate( fjac, njac, lhsa, lhsb, lhsc, rhsx )

      end module work_lhs
      

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine lhsinit(ni)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      use work_lhs
      implicit none

      integer ni

      integer i, m, n, jb

c---------------------------------------------------------------------
c     zero the whole left hand side for starters
c     set all diagonal values to 1. This is overkill, but convenient
c---------------------------------------------------------------------
      if (ni .gt. 0) goto 20

      do i = 0, 2, 2
!dir$ vector always
         do jb = 1, bsize
!dir$ unroll
            do m = 1, 5
!dir$ unroll
               do n = 1, 5
                  lhsa(jb,m,n,i) = 0.0d0
                  lhsb(jb,m,n,i) = 0.0d0
               end do
               lhsb(jb,m,m,i) = 1.0d0
            end do
         end do
      end do

      return

  20  continue
      do i = 0, ni, ni
!dir$ vector always
         do jb = 1, bsize
!dir$ unroll
            do m = 1, 5
!dir$ unroll
               do n = 1, 5
                  lhsc(jb,m,n,i) = 0.0d0
               end do
            end do
         end do
      end do

      return
      end

