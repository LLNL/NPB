c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  syncs module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module syncs

      use lu_data, only : isiz2

c---------------------------------------------------------------------
c  Flags used for thread synchronization for pipeline operation
c---------------------------------------------------------------------

      integer padim
      parameter (padim=16)
      integer isync(padim,0:isiz2), mthreadnum, iam
!$omp threadprivate( mthreadnum, iam )

      end module syncs


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine sync_init( jdim )

c---------------------------------------------------------------------
c   Initialize sync-related variables
c---------------------------------------------------------------------

      use syncs
      implicit none

      integer jdim

!$    integer, external :: omp_get_num_threads, omp_get_thread_num

      mthreadnum = 0
!$    mthreadnum = omp_get_num_threads() - 1
      if (mthreadnum .gt. jdim) mthreadnum = jdim
      iam = 0
!$    iam = omp_get_thread_num()
      if (iam .le. mthreadnum) isync(1,iam) = 0

      return
      end

c---------------------------------------------------------------------
c   Thread synchronization for pipeline operation
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine sync_left( ldmx, ldmy, ldmz, v )

c---------------------------------------------------------------------
c   Thread synchronization for pipeline operation
c---------------------------------------------------------------------

      use syncs
      implicit none

      integer ldmx, ldmy, ldmz
      double precision  v( 5, ldmx/2*2+1, ldmy/2*2+1, ldmz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer neigh, iv


      if (iam .gt. 0 .and. iam .le. mthreadnum) then
         neigh = iam - 1
!$omp atomic read
         iv = isync(1,neigh)
         do while (iv .eq. 0)
!$omp atomic read
            iv = isync(1,neigh)
         end do
!$omp atomic write
         isync(1,neigh) = 0
      endif
!$omp flush(isync,v)


      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine sync_right( ldmx, ldmy, ldmz, v )

c---------------------------------------------------------------------
c   Thread synchronization for pipeline operation
c---------------------------------------------------------------------

      use syncs
      implicit none

      integer ldmx, ldmy, ldmz
      double precision  v( 5, ldmx/2*2+1, ldmy/2*2+1, ldmz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer iv


!$omp flush(isync,v)
      if (iam .lt. mthreadnum) then
!$omp atomic read
         iv = isync(1,iam)
         do while (iv .eq. 1)
!$omp atomic read
            iv = isync(1,iam)
         end do
!$omp atomic write
         isync(1,iam) = 1
      endif


      return
      end
