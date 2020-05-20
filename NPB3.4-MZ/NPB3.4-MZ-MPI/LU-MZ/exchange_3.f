
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine exchange_3(g, nx, nxmax, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   compute the right hand side based on exact solution
c---------------------------------------------------------------------

      use lu_data
      use mpinpb

      implicit none

c---------------------------------------------------------------------
c  input parameters
c---------------------------------------------------------------------
      integer nx, nxmax, ny, nz
      double precision  g(5,-1:nxmax+2,ny,nz)

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer i, j, k
      integer ipos

      integer mid
      integer STATUS(MPI_STATUS_SIZE)


      if (timeron) call timer_start(t_exch)

c---------------------------------------------------------------------
c   communicate in the south and north directions
c---------------------------------------------------------------------
      if (north.ne.-1) then
          call MPI_IRECV( buf1,
     >                    10*ny*nz,
     >                    dp_type,
     >                    north,
     >                    from_n,
     >                    comm_setup,
     >                    mid,
     >                    IERROR )
      end if

c---------------------------------------------------------------------
c   send south
c---------------------------------------------------------------------
      if (south.ne.-1) then
          do k = 1,nz
            do j = 1,ny
              ipos = (k-1)*ny + j
              buf(1,1,ipos) = g(1,nx-1,j,k) 
              buf(2,1,ipos) = g(2,nx-1,j,k) 
              buf(3,1,ipos) = g(3,nx-1,j,k) 
              buf(4,1,ipos) = g(4,nx-1,j,k) 
              buf(5,1,ipos) = g(5,nx-1,j,k) 
              buf(1,2,ipos) = g(1,nx,j,k)
              buf(2,2,ipos) = g(2,nx,j,k)
              buf(3,2,ipos) = g(3,nx,j,k)
              buf(4,2,ipos) = g(4,nx,j,k)
              buf(5,2,ipos) = g(5,nx,j,k)
            end do
          end do

          call MPI_SEND( buf,
     >                   10*ny*nz,
     >                   dp_type,
     >                   south,
     >                   from_n,
     >                   comm_setup,
     >                   IERROR )
      end if

c---------------------------------------------------------------------
c   receive from north
c---------------------------------------------------------------------
      if (north.ne.-1) then
          call MPI_WAIT( mid, STATUS, IERROR )

          do k = 1,nz
            do j = 1,ny
              ipos = (k-1)*ny + j
              g(1,-1,j,k) = buf1(1,1,ipos)
              g(2,-1,j,k) = buf1(2,1,ipos)
              g(3,-1,j,k) = buf1(3,1,ipos)
              g(4,-1,j,k) = buf1(4,1,ipos)
              g(5,-1,j,k) = buf1(5,1,ipos)
              g(1, 0,j,k) = buf1(1,2,ipos)
              g(2, 0,j,k) = buf1(2,2,ipos)
              g(3, 0,j,k) = buf1(3,2,ipos)
              g(4, 0,j,k) = buf1(4,2,ipos)
              g(5, 0,j,k) = buf1(5,2,ipos)
            end do
          end do

      end if

      if (south.ne.-1) then
          call MPI_IRECV( buf1,
     >                    10*ny*nz,
     >                    dp_type,
     >                    south,
     >                    from_s,
     >                    comm_setup,
     >                    mid,
     >                    IERROR )
      end if

c---------------------------------------------------------------------
c   send north
c---------------------------------------------------------------------
      if (north.ne.-1) then
          do k = 1,nz
            do j = 1,ny
              ipos = (k-1)*ny + j
              buf(1,1,ipos) = g(1,1,j,k)
              buf(2,1,ipos) = g(2,1,j,k)
              buf(3,1,ipos) = g(3,1,j,k)
              buf(4,1,ipos) = g(4,1,j,k)
              buf(5,1,ipos) = g(5,1,j,k)
              buf(1,2,ipos) = g(1,2,j,k)
              buf(2,2,ipos) = g(2,2,j,k)
              buf(3,2,ipos) = g(3,2,j,k)
              buf(4,2,ipos) = g(4,2,j,k)
              buf(5,2,ipos) = g(5,2,j,k)
            end do
          end do

          call MPI_SEND( buf,
     >                   10*ny*nz,
     >                   dp_type,
     >                   north,
     >                   from_s,
     >                   comm_setup,
     >                   IERROR )
      end if

c---------------------------------------------------------------------
c   receive from south
c---------------------------------------------------------------------
      if (south.ne.-1) then
          call MPI_WAIT( mid, STATUS, IERROR )

          do k = 1,nz
            do j = 1,ny
              ipos = (k-1)*ny + j
              g(1,nx+1,j,k) = buf1(1,1,ipos)
              g(2,nx+1,j,k) = buf1(2,1,ipos)
              g(3,nx+1,j,k) = buf1(3,1,ipos)
              g(4,nx+1,j,k) = buf1(4,1,ipos)
              g(5,nx+1,j,k) = buf1(5,1,ipos)
              g(1,nx+2,j,k) = buf1(1,2,ipos)
              g(2,nx+2,j,k) = buf1(2,2,ipos)
              g(3,nx+2,j,k) = buf1(3,2,ipos)
              g(4,nx+2,j,k) = buf1(4,2,ipos)
              g(5,nx+2,j,k) = buf1(5,2,ipos)
            end do
          end do
      end if

      if (timeron) call timer_stop(t_exch)


      return
      end     
