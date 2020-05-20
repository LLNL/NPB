
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine exchange_6(g,jbeg,jfin1)

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
      integer jbeg, jfin1
      double precision  g(0:isiz2+1,0:isiz3+1)

c---------------------------------------------------------------------
c  local parameters
c---------------------------------------------------------------------
      integer k
      double precision  dum(isiz03)

      integer msgid3
      integer STATUS(MPI_STATUS_SIZE)
      integer IERROR



c---------------------------------------------------------------------
c   communicate in the east and west directions
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   receive from east
c---------------------------------------------------------------------
      if (jfin1.eq.ny) then
        call MPI_IRECV( dum,
     >                  nz,
     >                  dp_type,
     >                  east,
     >                  from_e,
     >                  MPI_COMM_WORLD,
     >                  msgid3,
     >                  IERROR )

        call MPI_WAIT( msgid3, STATUS, IERROR )

        do k = 1,nz
          g(ny+1,k) = dum(k)
        end do

      end if

c---------------------------------------------------------------------
c   send west
c---------------------------------------------------------------------
      if (jbeg.eq.1) then
        do k = 1,nz
          dum(k) = g(1,k)
        end do

        call MPI_SEND( dum,
     >                 nz,
     >                 dp_type,
     >                 west,
     >                 from_e,
     >                 MPI_COMM_WORLD,
     >                 IERROR )

      end if

      return
      end     
