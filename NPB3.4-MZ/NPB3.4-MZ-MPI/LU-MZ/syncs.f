c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  syncs module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module syncs

      use lu_data, only : problem_size

c---------------------------------------------------------------------
c  Flags used for thread synchronization for pipeline operation
c---------------------------------------------------------------------

      integer isync(0:problem_size)
      integer mthreadnum, iam
!$omp threadprivate( mthreadnum, iam )

      end module syncs


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine sync_init

c---------------------------------------------------------------------
c   Initialize synchronization variables
c---------------------------------------------------------------------

      use syncs
      implicit none

c---------------------------------------------------------------------
c---------------------------------------------------------------------

!$    integer, external :: omp_get_thread_num
!$    integer, external :: omp_get_num_threads

      mthreadnum = 0
!$    mthreadnum = omp_get_num_threads() - 1
      if (mthreadnum .gt. problem_size) mthreadnum = problem_size
      iam = 0
!$    iam = omp_get_thread_num()
      if (iam.le.mthreadnum) isync(iam) = 0

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine sync_left( nxmax, ny, nz, v )

c---------------------------------------------------------------------
c   Thread synchronization for pipeline operation
c---------------------------------------------------------------------

      use syncs
      implicit none

      integer nxmax, ny, nz
      double precision  v(5, -1:nxmax+2, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer neigh, iv


      if (iam .gt. 0 .and. iam .le. mthreadnum) then
         neigh = iam - 1
!$omp atomic read
         iv = isync(neigh)
         do while (iv .eq. 0)
!$omp atomic read
            iv = isync(neigh)
         end do
!$omp atomic write
         isync(neigh) = 0
      endif
!$omp flush(isync,v)


      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine sync_right( nxmax, ny, nz, v )

c---------------------------------------------------------------------
c   Thread synchronization for pipeline operation
c---------------------------------------------------------------------

      use syncs
      implicit none

      integer nxmax, ny, nz
      double precision  v(5, -1:nxmax+2, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer iv


!$omp flush(isync,v)
      if (iam .lt. mthreadnum) then
!$omp atomic read
         iv = isync(iam)
         do while (iv .eq. 1)
!$omp atomic read
            iv = isync(iam)
         end do
!$omp atomic write
         isync(iam) = 1
      endif


      return
      end
