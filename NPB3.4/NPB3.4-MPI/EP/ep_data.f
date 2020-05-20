c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  ep_data module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------
 
      module ep_data

c---------------------------------------------------------------------
c  The following include file is generated automatically by the
c  "setparams" utility, which defines the problem size 'm'
c---------------------------------------------------------------------

      include 'npbparams.h'

c---------------------------------------------------------------------
c   M is the Log_2 of the number of complex pairs of uniform (0, 1) random
c   numbers.  MK is the Log_2 of the size of each batch of uniform random
c   numbers.  MK can be set for convenience on a given system, since it does
c   not affect the results.
c---------------------------------------------------------------------
      integer    mk, mm, nn, nk, nq
      parameter (mk = 16, mm = m - mk, nn = 2 ** mm,
     >           nk = 2 ** mk, nq = 10)

      double precision a, s
      parameter (a = 1220703125.d0, s = 271828183.d0)

c ... storage
      double precision x(2*nk), q(0:nq-1), qq(10000)

c ... timer constants
      integer    t_total, t_gpairs, t_randn, t_rcomm, t_last
      parameter (t_total=1, t_gpairs=2, t_randn=3, t_rcomm=4, t_last=4)

      end module ep_data

