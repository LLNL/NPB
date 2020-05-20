c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  ft_data module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module ft_data

      include 'npbparams.h'

c total number of grid points in floating point number
      double precision ntotal_f
      parameter (ntotal_f=dble(nx)*ny*nz)

c total dimension scaled by the number of processes
      integer ntdivnp


      double precision seed, a, pi, alpha
      parameter (seed = 314159265.d0, a = 1220703125.d0, 
     >  pi = 3.141592653589793238d0, alpha=1.0d-6)

c roots of unity array
c relies on x being largest dimension?
      double complex, allocatable :: u(:)


c for checksum data
      double complex sums(0:niter_default)

c number of iterations
      integer niter

c other stuff
      logical debug, debugsynch


c--------------------------------------------------------------------
c Cache blocking params. These values are good for most
c RISC processors.  
c FFT parameters:
c  fftblock controls how many ffts are done at a time. 
c  The default is appropriate for most cache-based machines
c  On vector machines, the FFT can be vectorized with vector
c  length equal to the block size, so the block size should
c  be as large as possible. This is the size of the smallest
c  dimension of the problem: 128 for class A, 256 for class B and
c  512 for class C.
c Transpose parameters:
c  transblock is the blocking factor for the transposes when there
c  is a 1-D layout. On vector machines it should probably be
c  large (largest dimension of the problem).
c--------------------------------------------------------------------

      integer fftblock_default, fftblockpad_default
      parameter (fftblock_default=16, fftblockpad_default=18)
      integer transblock, transblockpad
      parameter(transblock=32, transblockpad=34)
      
      integer fftblock, fftblockpad


c--------------------------------------------------------------------
c 2D processor array -> 2D grid decomposition (by pencils)
c If processor array is 1xN or -> 1D grid decomposition (by planes)
c If processor array is 1x1 -> 0D grid decomposition
c For simplicity, do not treat Nx1 (np2 = 1) specially
c--------------------------------------------------------------------
      integer np1, np2

c basic decomposition strategy
      integer layout_type
      integer layout_0D, layout_1D, layout_2D
      parameter (layout_0D = 0, layout_1D = 1, layout_2D = 2)

c--------------------------------------------------------------------
c There are basically three stages
c 1: x-y-z layout
c 2: after x-transform (before y)
c 3: after y-transform (before z)
c The computation proceeds logically as

c set up initial conditions
c fftx(1)
c transpose (1->2)
c ffty(2)
c transpose (2->3)
c fftz(3)
c time evolution
c fftz(3)
c transpose (3->2)
c ffty(2)
c transpose (2->1)
c fftx(1)
c compute residual(1)

c for the 0D, 1D, 2D strategies, the layouts look like xxx
c        
c            0D        1D        2D
c 1:        xyz       xyz       xyz
c 2:        xyz       xyz       yxz
c 3:        xyz       zyx       zxy
c--------------------------------------------------------------------

c the array dimensions are stored in dims(coord, phase)
      integer dims(3, 3)
      integer xstart(3), ystart(3), zstart(3)
      integer xend(3), yend(3), zend(3)

c--------------------------------------------------------------------
c Timing constants
c--------------------------------------------------------------------
      integer T_total, T_setup, T_fft, T_evolve, T_checksum, 
     >        T_fftlow, T_fftcopy, T_transpose, 
     >        T_transxzloc, T_transxzglo, T_transxzfin, 
     >        T_transxyloc, T_transxyglo, T_transxyfin, 
     >        T_synch, T_init, T_max
      parameter (T_total = 1, T_setup = 2, T_fft = 3, 
     >           T_evolve = 4, T_checksum = 5, 
     >           T_fftlow = 6, T_fftcopy = 7, T_transpose = 8,
     >           T_transxzloc = 9, T_transxzglo = 10, T_transxzfin = 11, 
     >           T_transxyloc = 12, T_transxyglo = 13, 
     >           T_transxyfin = 14,  T_synch = 15, T_init = 16,
     >           T_max = 16)

      logical timers_enabled

c--------------------------------------------------------------------
c external functions
c--------------------------------------------------------------------
      double precision, external :: randlc, timer_read
      integer, external ::          ilog2

      end module ft_data


c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  ft_fields module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module ft_fields

c---------------------------------------------------------------------
c u0, u1, u2 are the main arrays in the problem. 
c Depending on the decomposition, these arrays will have different 
c dimensions. To accomodate all possibilities, we allocate them as 
c one-dimensional arrays and pass them to subroutines for different 
c views
c  - u0 contains the initial (transformed) initial condition
c  - u1 and u2 are working arrays
c---------------------------------------------------------------------
      double complex, allocatable ::
     >                 u0(:), u1(:), u2(:)
      double precision, allocatable ::
     >                 twiddle(:)

      end module ft_fields


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine alloc_space

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c allocate space dynamically for data arrays
c---------------------------------------------------------------------

      use ft_data
      use ft_fields
      use mpinpb

      implicit none

      integer ios, ierr


      ntdivnp = ((nx*ny)/np_min)*nz

c---------------------------------------------------------------------
c Padding+3 is to avoid accidental cache problems, 
c since all array sizes are powers of two.
c---------------------------------------------------------------------
      allocate (
     >          u0     (ntdivnp+3), 
     >          u1     (ntdivnp+3), 
     >          u2     (ntdivnp+3),
     >          twiddle(ntdivnp),
     >          u      (maxdim),
     >          stat = ios)

      if (ios .ne. 0) then
         write(*,*) 'Error encountered in allocating space'
         call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
         stop
      endif

      return
      end

