c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  ft_data module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module ft_data

      include 'npbparams.h'

c total number of grid points with padding
      integer(kind2) nxp, ntotalp
      parameter (nxp=nx+1)
      parameter (ntotalp=nxp*ny*nz)
      double precision ntotal_f
      parameter (ntotal_f=dble(nx)*ny*nz)


c If processor array is 1x1 -> 0D grid decomposition


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

      integer fftblock_default, fftblockpad_default
c      parameter (fftblock_default=16, fftblockpad_default=18)
      parameter (fftblock_default=32, fftblockpad_default=33)
      
      integer fftblock, fftblockpad

c we need a bunch of logic to keep track of how
c arrays are laid out. 


c Note: this serial version is the derived from the parallel 0D case
c of the ft NPB.
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

c the array dimensions are stored in dims(coord, phase)
      integer dims(3)

      integer T_total, T_setup, T_fft, T_evolve, T_checksum, 
     >        T_fftx, T_ffty,
     >        T_fftz, T_max
      parameter (T_total = 1, T_setup = 2, T_fft = 3, 
     >           T_evolve = 4, T_checksum = 5, 
     >           T_fftx = 6,
     >           T_ffty = 7,
     >           T_fftz = 8, T_max = 8)



      logical timers_enabled


      external timer_read
      double precision timer_read
      external ilog2
      integer ilog2

      external randlc
      double precision randlc


c other stuff
      logical debug, debugsynch

      double precision seed, a, pi, alpha
      parameter (seed = 314159265.d0, a = 1220703125.d0, 
     >  pi = 3.141592653589793238d0, alpha=1.0d-6)


c roots of unity array
c relies on x being largest dimension?
      double complex u(nxp)


c for checksum data
      double complex sums(0:niter_default)

c number of iterations
      integer niter


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
c  - twiddle contains exponents for the time evolution operator. 
c---------------------------------------------------------------------

      double complex, allocatable ::
     >                 u0(:), pad1(:),
     >                 u1(:), pad2(:)
c     >                 u2(:)
      double precision, allocatable :: twiddle(:)
c---------------------------------------------------------------------
c Large arrays are in module so that they are allocated on the
c heap rather than the stack. This module is not
c referenced directly anywhere else. Padding is to avoid accidental 
c cache problems, since all array sizes are powers of two.
c---------------------------------------------------------------------


      end module ft_fields


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine alloc_space

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c allocate space dynamically for data arrays
c---------------------------------------------------------------------

      use ft_data, only : ntotalp
      use ft_fields

      implicit none

      integer ios


      allocate ( 
     >          u0(ntotalp), pad1(3),
     >          u1(ntotalp), pad2(3),
c     >          u2(ntotalp),
     >          twiddle(ntotalp),
     >          stat = ios)

      if (ios .ne. 0) then
         write(*,*) 'Error encountered in allocating space'
         stop
      endif

      return
      end


