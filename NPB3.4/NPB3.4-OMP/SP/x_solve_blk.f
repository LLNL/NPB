
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       subroutine x_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c this function performs the solution of the approximate factorization
c step in the x-direction for all five matrix components
c simultaneously. The Thomas algorithm is employed to solve the
c systems for the x-lines. Boundary conditions are non-periodic
c---------------------------------------------------------------------

       use sp_data
       use work_lhs

       implicit none

       integer i, j, k, i1, i2, jj, jb, jm
       double precision  ru1, fac1, fac2


c---------------------------------------------------------------------
c---------------------------------------------------------------------

       if (timeron) call timer_start(t_xsolve)
!$omp parallel default(shared) private(i,j,k,i1,i2,jj,jb,jm,
!$omp&    ru1,fac1,fac2)

       call lhsinit(nx2+1)

!$omp do collapse(2)
       do  k = 1, nz2
       do  jj = 1, ny2, bsize
          jm = min(bsize, ny2 - jj + 1)

c---------------------------------------------------------------------
c To improve cache utilization, copy a slab of rhs to temp array  
c---------------------------------------------------------------------
          do  i = 0, grid_points(1)-1
             do  jb = 1, bsize
                j = min(jb,jm) + jj - 1
                rhsx(jb,1,i) = rhs(1,i,j,k)
                rhsx(jb,2,i) = rhs(2,i,j,k)
                rhsx(jb,3,i) = rhs(3,i,j,k)
                rhsx(jb,4,i) = rhs(4,i,j,k)
                rhsx(jb,5,i) = rhs(5,i,j,k)
             end do
          end do

c---------------------------------------------------------------------
c Computes the left hand side for the three x-factors  
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      first fill the lhs for the u-eigenvalue                   
c---------------------------------------------------------------------
          do  i = 0, grid_points(1)-1
             do  jb = 1, bsize
                j = min(jb,jm) + jj - 1
                ru1 = c3c4*rho_i(i,j,k)
                cv(jb,i) = us(i,j,k)
                rhov(jb,i) = dmax1(dx2+con43*ru1, 
     >                          dx5+c1c5*ru1,
     >                          dxmax+ru1,
     >                          dx1)
             end do
          end do

          do  i = 1, nx2
             do  jb = 1, bsize
                lhs(jb,1,i) =  0.0d0
                lhs(jb,2,i) = -dttx2 * cv(jb,i-1) - dttx1 * rhov(jb,i-1)
                lhs(jb,3,i) =  1.0d0 + c2dttx1 * rhov(jb,i)
                lhs(jb,4,i) =  dttx2 * cv(jb,i+1) - dttx1 * rhov(jb,i+1)
                lhs(jb,5,i) =  0.0d0
             end do
          end do

c---------------------------------------------------------------------
c      add fourth order dissipation                             
c---------------------------------------------------------------------

          do  jb = 1, bsize
             i = 1
             lhs(jb,3,i) = lhs(jb,3,i) + comz5
             lhs(jb,4,i) = lhs(jb,4,i) - comz4
             lhs(jb,5,i) = lhs(jb,5,i) + comz1
  
             i = 2
             lhs(jb,2,i) = lhs(jb,2,i) - comz4
             lhs(jb,3,i) = lhs(jb,3,i) + comz6
             lhs(jb,4,i) = lhs(jb,4,i) - comz4
             lhs(jb,5,i) = lhs(jb,5,i) + comz1
          end do

          do   i=3, grid_points(1)-4
             do  jb = 1, bsize
                lhs(jb,1,i) = lhs(jb,1,i) + comz1
                lhs(jb,2,i) = lhs(jb,2,i) - comz4
                lhs(jb,3,i) = lhs(jb,3,i) + comz6
                lhs(jb,4,i) = lhs(jb,4,i) - comz4
                lhs(jb,5,i) = lhs(jb,5,i) + comz1
             end do
          end do

          do  jb = 1, bsize
             i = grid_points(1)-3
             lhs(jb,1,i) = lhs(jb,1,i) + comz1
             lhs(jb,2,i) = lhs(jb,2,i) - comz4
             lhs(jb,3,i) = lhs(jb,3,i) + comz6
             lhs(jb,4,i) = lhs(jb,4,i) - comz4

             i = grid_points(1)-2
             lhs(jb,1,i) = lhs(jb,1,i) + comz1
             lhs(jb,2,i) = lhs(jb,2,i) - comz4
             lhs(jb,3,i) = lhs(jb,3,i) + comz5
          end do

c---------------------------------------------------------------------
c      subsequently, fill the other factors (u+c), (u-c) by adding to 
c      the first  
c---------------------------------------------------------------------
          do   i = 1, nx2
             do  jb = 1, bsize
                j = min(jb,jm) + jj - 1
                lhsp(jb,1,i) = lhs(jb,1,i)
                lhsp(jb,2,i) = lhs(jb,2,i) - 
     >                            dttx2 * speed(i-1,j,k)
                lhsp(jb,3,i) = lhs(jb,3,i)
                lhsp(jb,4,i) = lhs(jb,4,i) + 
     >                            dttx2 * speed(i+1,j,k)
                lhsp(jb,5,i) = lhs(jb,5,i)
                lhsm(jb,1,i) = lhs(jb,1,i)
                lhsm(jb,2,i) = lhs(jb,2,i) + 
     >                            dttx2 * speed(i-1,j,k)
                lhsm(jb,3,i) = lhs(jb,3,i)
                lhsm(jb,4,i) = lhs(jb,4,i) - 
     >                            dttx2 * speed(i+1,j,k)
                lhsm(jb,5,i) = lhs(jb,5,i)
             end do
          end do

c---------------------------------------------------------------------
c                          FORWARD ELIMINATION  
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      perform the Thomas algorithm; first, FORWARD ELIMINATION     
c---------------------------------------------------------------------

          do    i = 0, grid_points(1)-3
             i1 = i  + 1
             i2 = i  + 2
             do  jb = 1, bsize
                fac1      = 1.d0/lhs(jb,3,i)
                lhs(jb,4,i)  = fac1*lhs(jb,4,i)
                lhs(jb,5,i)  = fac1*lhs(jb,5,i)
                rhsx(jb,1,i) = fac1*rhsx(jb,1,i)
                rhsx(jb,2,i) = fac1*rhsx(jb,2,i)
                rhsx(jb,3,i) = fac1*rhsx(jb,3,i)
                lhs(jb,3,i1) = lhs(jb,3,i1) -
     >                         lhs(jb,2,i1)*lhs(jb,4,i)
                lhs(jb,4,i1) = lhs(jb,4,i1) -
     >                         lhs(jb,2,i1)*lhs(jb,5,i)
                rhsx(jb,1,i1) = rhsx(jb,1,i1) -
     >                         lhs(jb,2,i1)*rhsx(jb,1,i)
                rhsx(jb,2,i1) = rhsx(jb,2,i1) -
     >                         lhs(jb,2,i1)*rhsx(jb,2,i)
                rhsx(jb,3,i1) = rhsx(jb,3,i1) -
     >                         lhs(jb,2,i1)*rhsx(jb,3,i)
                lhs(jb,2,i2) = lhs(jb,2,i2) -
     >                         lhs(jb,1,i2)*lhs(jb,4,i)
                lhs(jb,3,i2) = lhs(jb,3,i2) -
     >                         lhs(jb,1,i2)*lhs(jb,5,i)
                rhsx(jb,1,i2) = rhsx(jb,1,i2) -
     >                         lhs(jb,1,i2)*rhsx(jb,1,i)
                rhsx(jb,2,i2) = rhsx(jb,2,i2) -
     >                         lhs(jb,1,i2)*rhsx(jb,2,i)
                rhsx(jb,3,i2) = rhsx(jb,3,i2) -
     >                         lhs(jb,1,i2)*rhsx(jb,3,i)
             end do
          end do

c---------------------------------------------------------------------
c      The last two rows in this grid block are a bit different, 
c      since they do not have two more rows available for the
c      elimination of off-diagonal entries
c---------------------------------------------------------------------

          i  = grid_points(1)-2
          i1 = grid_points(1)-1
          do  jb = 1, bsize
             fac1      = 1.d0/lhs(jb,3,i)
             lhs(jb,4,i)  = fac1*lhs(jb,4,i)
             lhs(jb,5,i)  = fac1*lhs(jb,5,i)
             rhsx(jb,1,i) = fac1*rhsx(jb,1,i)
             rhsx(jb,2,i) = fac1*rhsx(jb,2,i)
             rhsx(jb,3,i) = fac1*rhsx(jb,3,i)
             lhs(jb,3,i1) = lhs(jb,3,i1) -
     >                      lhs(jb,2,i1)*lhs(jb,4,i)
             lhs(jb,4,i1) = lhs(jb,4,i1) -
     >                      lhs(jb,2,i1)*lhs(jb,5,i)
             rhsx(jb,1,i1) = rhsx(jb,1,i1) -
     >                      lhs(jb,2,i1)*rhsx(jb,1,i)
             rhsx(jb,2,i1) = rhsx(jb,2,i1) -
     >                      lhs(jb,2,i1)*rhsx(jb,2,i)
             rhsx(jb,3,i1) = rhsx(jb,3,i1) -
     >                      lhs(jb,2,i1)*rhsx(jb,3,i)
c---------------------------------------------------------------------
c            scale the last row immediately 
c---------------------------------------------------------------------
             fac2             = 1.d0/lhs(jb,3,i1)
             rhsx(jb,1,i1) = fac2*rhsx(jb,1,i1)
             rhsx(jb,2,i1) = fac2*rhsx(jb,2,i1)
             rhsx(jb,3,i1) = fac2*rhsx(jb,3,i1)
          end do

c---------------------------------------------------------------------
c      do the u+c and the u-c factors                 
c---------------------------------------------------------------------

          do    i = 0, grid_points(1)-3
             i1 = i  + 1
             i2 = i  + 2
             do  jb = 1, bsize
                fac1       = 1.d0/lhsp(jb,3,i)
                lhsp(jb,4,i)  = fac1*lhsp(jb,4,i)
                lhsp(jb,5,i)  = fac1*lhsp(jb,5,i)
                rhsx(jb,4,i)  = fac1*rhsx(jb,4,i)
                lhsp(jb,3,i1) = lhsp(jb,3,i1) -
     >                        lhsp(jb,2,i1)*lhsp(jb,4,i)
                lhsp(jb,4,i1) = lhsp(jb,4,i1) -
     >                        lhsp(jb,2,i1)*lhsp(jb,5,i)
                rhsx(jb,4,i1) = rhsx(jb,4,i1) -
     >                        lhsp(jb,2,i1)*rhsx(jb,4,i)
                lhsp(jb,2,i2) = lhsp(jb,2,i2) -
     >                        lhsp(jb,1,i2)*lhsp(jb,4,i)
                lhsp(jb,3,i2) = lhsp(jb,3,i2) -
     >                        lhsp(jb,1,i2)*lhsp(jb,5,i)
                rhsx(jb,4,i2) = rhsx(jb,4,i2) -
     >                        lhsp(jb,1,i2)*rhsx(jb,4,i)
                fac1       = 1.d0/lhsm(jb,3,i)
                lhsm(jb,4,i)  = fac1*lhsm(jb,4,i)
                lhsm(jb,5,i)  = fac1*lhsm(jb,5,i)
                rhsx(jb,5,i)  = fac1*rhsx(jb,5,i)
                lhsm(jb,3,i1) = lhsm(jb,3,i1) -
     >                        lhsm(jb,2,i1)*lhsm(jb,4,i)
                lhsm(jb,4,i1) = lhsm(jb,4,i1) -
     >                        lhsm(jb,2,i1)*lhsm(jb,5,i)
                rhsx(jb,5,i1) = rhsx(jb,5,i1) -
     >                        lhsm(jb,2,i1)*rhsx(jb,5,i)
                lhsm(jb,2,i2) = lhsm(jb,2,i2) -
     >                        lhsm(jb,1,i2)*lhsm(jb,4,i)
                lhsm(jb,3,i2) = lhsm(jb,3,i2) -
     >                        lhsm(jb,1,i2)*lhsm(jb,5,i)
                rhsx(jb,5,i2) = rhsx(jb,5,i2) -
     >                        lhsm(jb,1,i2)*rhsx(jb,5,i)
             end do
          end do

c---------------------------------------------------------------------
c         And again the last two rows separately
c---------------------------------------------------------------------
          i  = grid_points(1)-2
          i1 = grid_points(1)-1
          do  jb = 1, bsize
             fac1       = 1.d0/lhsp(jb,3,i)
             lhsp(jb,4,i)  = fac1*lhsp(jb,4,i)
             lhsp(jb,5,i)  = fac1*lhsp(jb,5,i)
             rhsx(jb,4,i)  = fac1*rhsx(jb,4,i)
             lhsp(jb,3,i1) = lhsp(jb,3,i1) -
     >                      lhsp(jb,2,i1)*lhsp(jb,4,i)
             lhsp(jb,4,i1) = lhsp(jb,4,i1) -
     >                      lhsp(jb,2,i1)*lhsp(jb,5,i)
             rhsx(jb,4,i1) = rhsx(jb,4,i1) -
     >                      lhsp(jb,2,i1)*rhsx(jb,4,i)
             fac1       = 1.d0/lhsm(jb,3,i)
             lhsm(jb,4,i)  = fac1*lhsm(jb,4,i)
             lhsm(jb,5,i)  = fac1*lhsm(jb,5,i)
             rhsx(jb,5,i)  = fac1*rhsx(jb,5,i)
             lhsm(jb,3,i1) = lhsm(jb,3,i1) -
     >                      lhsm(jb,2,i1)*lhsm(jb,4,i)
             lhsm(jb,4,i1) = lhsm(jb,4,i1) -
     >                      lhsm(jb,2,i1)*lhsm(jb,5,i)
             rhsx(jb,5,i1) = rhsx(jb,5,i1) -
     >                      lhsm(jb,2,i1)*rhsx(jb,5,i)
c---------------------------------------------------------------------
c               Scale the last row immediately
c---------------------------------------------------------------------
             rhsx(jb,4,i1) = rhsx(jb,4,i1)/lhsp(jb,3,i1)
             rhsx(jb,5,i1) = rhsx(jb,5,i1)/lhsm(jb,3,i1)
          end do


c---------------------------------------------------------------------
c                         BACKSUBSTITUTION 
c---------------------------------------------------------------------


          i  = grid_points(1)-2
          i1 = grid_points(1)-1
          do  jb = 1, bsize
             rhsx(jb,1,i) = rhsx(jb,1,i) -
     >                             lhs(jb,4,i)*rhsx(jb,1,i1)
             rhsx(jb,2,i) = rhsx(jb,2,i) -
     >                             lhs(jb,4,i)*rhsx(jb,2,i1)
             rhsx(jb,3,i) = rhsx(jb,3,i) -
     >                             lhs(jb,4,i)*rhsx(jb,3,i1)

             rhsx(jb,4,i) = rhsx(jb,4,i) -
     >                          lhsp(jb,4,i)*rhsx(jb,4,i1)
             rhsx(jb,5,i) = rhsx(jb,5,i) -
     >                          lhsm(jb,4,i)*rhsx(jb,5,i1)
          end do

c---------------------------------------------------------------------
c      The first three factors
c---------------------------------------------------------------------
          do    i = grid_points(1)-3, 0, -1
             i1 = i  + 1
             i2 = i  + 2
             do  jb = 1, bsize
                rhsx(jb,1,i) = rhsx(jb,1,i) - 
     >                          lhs(jb,4,i)*rhsx(jb,1,i1) -
     >                          lhs(jb,5,i)*rhsx(jb,1,i2)
                rhsx(jb,2,i) = rhsx(jb,2,i) - 
     >                          lhs(jb,4,i)*rhsx(jb,2,i1) -
     >                          lhs(jb,5,i)*rhsx(jb,2,i2)
                rhsx(jb,3,i) = rhsx(jb,3,i) - 
     >                          lhs(jb,4,i)*rhsx(jb,3,i1) -
     >                          lhs(jb,5,i)*rhsx(jb,3,i2)

c---------------------------------------------------------------------
c      And the remaining two
c---------------------------------------------------------------------
                rhsx(jb,4,i) = rhsx(jb,4,i) - 
     >                          lhsp(jb,4,i)*rhsx(jb,4,i1) -
     >                          lhsp(jb,5,i)*rhsx(jb,4,i2)
                rhsx(jb,5,i) = rhsx(jb,5,i) - 
     >                          lhsm(jb,4,i)*rhsx(jb,5,i1) -
     >                          lhsm(jb,5,i)*rhsx(jb,5,i2)
             end do
          end do

          do  jb = 1, jm
             j = jb + jj - 1
             do  i = 0, grid_points(1)-1
                rhs(1,i,j,k) = rhsx(jb,1,i)
                rhs(2,i,j,k) = rhsx(jb,2,i)
                rhs(3,i,j,k) = rhsx(jb,3,i)
                rhs(4,i,j,k) = rhsx(jb,4,i)
                rhs(5,i,j,k) = rhsx(jb,5,i)
             end do
          end do

       end do
       end do
!$omp end do nowait
!$omp end parallel
       if (timeron) call timer_stop(t_xsolve)

c---------------------------------------------------------------------
c      Do the block-diagonal inversion          
c---------------------------------------------------------------------
       call ninvr

       return
       end


