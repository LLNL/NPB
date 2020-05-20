c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  mg_data module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module mg_data

c---------------------------------------------------------------------
c  Parameter lm is the log-base2 of the edge size max for
c  the partition on a given node, so must be changed either
c  to save space (if running a small case) or made bigger for larger 
c  cases, for example, 512^3. Thus lm=7 means that the largest dimension 
c  of a partition that can be solved on a node is 2^7 = 128. lm is set 
c  automatically in npbparams.h
c  Parameters ndim1, ndim2, ndim3 are the local problem dimensions. 
c---------------------------------------------------------------------

      include 'npbparams.h'

      ! partitioned size in each dimension
      integer ndim1, ndim2, ndim3

      ! log of maximum dimension on a node
      integer lm

      integer nm      ! actual dimension including ghost cells for communications
     >      , nv      ! size of rhs array
     >      , nr      ! size of residual array
     >      , nm2     ! size of communication buffer
     >      , maxlevel! maximum number of levels
      parameter (maxlevel = lt_default+1)


      integer nprocs_total

      integer maxprocs
      parameter( maxprocs = 131072 )  ! this is the upper proc limit that 
                                      ! the current "nr" parameter can handle
c---------------------------------------------------------------------
      integer nbr(3,-1:1,maxlevel), msg_type(3,-1:1)
      integer  msg_id(3,-1:1,2),nx(maxlevel),ny(maxlevel),nz(maxlevel)

      character class

      integer debug_vec(0:7)

      integer ir(maxlevel), m1(maxlevel), m2(maxlevel), m3(maxlevel)
      integer lt, lb

      logical dead(maxlevel), give_ex(3,maxlevel), take_ex(3,maxlevel)

c ... grid
      integer  is1, is2, is3, ie1, ie2, ie3

c---------------------------------------------------------------------
c  Set at m=1024, can handle cases up to 1024^3 case
c---------------------------------------------------------------------
      integer m
c      parameter( m=1037 )

      double precision, allocatable ::
     >        buff(:,:)

c---------------------------------------------------------------------
c  Timing constants
c---------------------------------------------------------------------
      integer t_bench, t_init, t_psinv, t_resid, t_rprj3, t_interp, 
     >        t_norm2u3, t_comm3, t_rcomm, t_last
      parameter (t_bench=1, t_init=2, t_psinv=3, t_resid=4, t_rprj3=5,  
     >        t_interp=6, t_norm2u3=7, t_comm3=8, 
     >        t_rcomm=9, t_last=9)

      logical timeron


      end module mg_data


c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  mg_fields module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module mg_fields

c---------------------------------------------------------------------------c
c These are major data arrays and can be quite large.
c They are always passed as subroutine args.
c---------------------------------------------------------------------------c
      double precision, allocatable :: u(:), v(:), r(:)

      double precision  a(0:3),c(0:3)

      end module mg_fields


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine alloc_space

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c allocate space dynamically for data arrays
c---------------------------------------------------------------------

      use mg_data
      use mg_fields
      use mpinpb

      implicit none

      integer ios, ierr
      integer log2_size, log_p


c---------------------------------------------------------------------
c set up dimension parameters after partition
c---------------------------------------------------------------------
      log_p  = log(float(nprocs)+0.0001)/log(2.0)

      ! lt is log of largest total dimension
      log2_size = lt_default

      ! log of maximum dimension on a node
      lm = log2_size - log_p/3
      ndim1 = lm
      ndim3 = log2_size - (log_p+2)/3
      ndim2 = log2_size - (log_p+1)/3

      ! array size parameters
      nm = 2+2**lm
      nv = (2+2**ndim1)*(2+2**ndim2)*(2+2**ndim3)
      nm2= 2*nm*nm
      nr = (8*(nv+nm**2+5*nm+14*lt_default-7*lm))/7
      m  = nm + 1

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      allocate (
     >          u(nr),
     >          v(nv),
     >          r(nr),
     >          buff(nm2,4),
     >          stat = ios)

      if (ios .ne. 0) then
         write(*,*) 'Error encountered in allocating space'
         call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
         stop
      endif

      return
      end

