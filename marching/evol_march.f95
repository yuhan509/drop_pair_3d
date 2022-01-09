SUBROUTINE expl_fwd_elr(ntheta, nphi, dt, ugrids, xgrid, ygrid, zgrid) 
   implicit real *8 (a-h,o-z)
   integer :: ntheta, nphi
   real *8 :: dt, tol
   real *8:: xgrid(nphi,ntheta),ygrid(nphi,ntheta),zgrid(nphi,ntheta)
   real *8:: ugrids(nphi,ntheta,3)

   real *8, ALLOCATABLE :: xt(:,:), yt(:,:), zt(:,:)
   allocate(xt(nphi,ntheta),yt(nphi,ntheta),zt(nphi,ntheta))
   
   do j = 1,nphi
     do i = 1,ntheta
        xgrid(j,i) = dt * ugrids(j,i,1) + xgrid(j,i)
        ygrid(j,i) = dt * ugrids(j,i,2) + ygrid(j,i)
        zgrid(j,i) = dt * ugrids(j,i,3) + zgrid(j,i)
     enddo               
  enddo
   
   END SUBROUTINE


SUBROUTINE expl_midp_pair(istep, dt,  &
u1grids, x1grid, y1grid, z1grid,      &
u2grids, x2grid, y2grid, z2grid) 

use mod_orgHelper
use mod_assign

integer,intent(IN) :: istep
real *8,intent(INOUT) :: dt
real *8,intent(INOUT) :: x1grid(nphi,ntheta),y1grid(nphi,ntheta),z1grid(nphi,ntheta)
real *8,intent(INOUT) :: x2grid(nphi,ntheta),y2grid(nphi,ntheta),z2grid(nphi,ntheta)
real *8,intent(IN) :: u1grids(nphi,ntheta,3),u2grids(nphi,ntheta,3)

integer :: dtShrinked, num_run
real *8, ALLOCATABLE, dimension(:,:) :: x1h, y1h, z1h, x1e, y1e, z1e, &
x1m, y1m, z1m, x2h, y2h, z2h, x2e, y2e, z2e,x2m, y2m, z2m   
real *8, ALLOCATABLE, dimension(:,:,:) :: u1h, u2h
complex *16, ALLOCATABLE, dimension(:,:) :: x1mpole, y1mpole, z1mpole, &
x2mpole, y2mpole, z2mpole
real *8, ALLOCATABLE, dimension(:,:) :: h1grid,ws1grid,h2grid,ws2grid
real *8, ALLOCATABLE, dimension(:,:,:) :: n1grid, n2grid
complex *16, ALLOCATABLE, dimension(:,:) :: h1mpole, ws1mpole, h2mpole, ws2mpole
complex *16, ALLOCATABLE, dimension(:,:,:) :: n1mpole, n2mpole
complex *16, ALLOCATABLE, dimension(:,:) :: en1mpole, en2mpole
real*8 :: max1, max2, newDt1, newDt2

allocate(x1h(nphi,ntheta),y1h(nphi,ntheta),z1h(nphi,ntheta))
allocate(x1e(nphi,ntheta),y1e(nphi,ntheta),z1e(nphi,ntheta))
allocate(x1m(nphi,ntheta),y1m(nphi,ntheta),z1m(nphi,ntheta))
allocate(x2h(nphi,ntheta),y2h(nphi,ntheta),z2h(nphi,ntheta))
allocate(x2e(nphi,ntheta),y2e(nphi,ntheta),z2e(nphi,ntheta))
allocate(x2m(nphi,ntheta),y2m(nphi,ntheta),z2m(nphi,ntheta))
allocate(u1h(nphi,ntheta,3),u2h(nphi,ntheta,3))

allocate(x1mpole(0:nterms,0:nterms),y1mpole(0:nterms,0:nterms),z1mpole(0:nterms,0:nterms))
allocate(x2mpole(0:nterms,0:nterms),y2mpole(0:nterms,0:nterms),z2mpole(0:nterms,0:nterms))

allocate(h1grid(nphi,ntheta),ws1grid(nphi,ntheta),h2grid(nphi,ntheta),ws2grid(nphi,ntheta))
allocate(n1grid(nphi,ntheta,3), n2grid(nphi,ntheta,3))
allocate(h1mpole(0:nterms,0:nterms), ws1mpole(0:nterms,0:nterms),  &
h2mpole(0:nterms,0:nterms), ws2mpole(0:nterms,0:nterms))
allocate(n1mpole(0:nterms,0:nterms,3), n2mpole(0:nterms,0:nterms,3))
allocate(en1mpole(0:nterms,0:nterms),en2mpole(0:nterms,0:nterms))

num_run = 0
dtShrinked = 1
do while (dtShrinked == 1)
  !!! prepare for error estimation
  do j = 1,nphi
     do i = 1,ntheta
        x1e(j,i) = dt * u1grids(j,i,1) + x1grid(j,i)
        y1e(j,i) = dt * u1grids(j,i,2) + y1grid(j,i)
        z1e(j,i) = dt * u1grids(j,i,3) + z1grid(j,i)
        x2e(j,i) = dt * u2grids(j,i,1) + x2grid(j,i)
        y2e(j,i) = dt * u2grids(j,i,2) + y2grid(j,i)
        z2e(j,i) = dt * u2grids(j,i,3) + z2grid(j,i) 
     enddo
  enddo
  !! marching to half step position
  do j = 1,nphi
     do i = 1,ntheta
        x1h(j,i) = dt/2 * u1grids(j,i,1) + x1grid(j,i)
        y1h(j,i) = dt/2 * u1grids(j,i,2) + y1grid(j,i)
        z1h(j,i) = dt/2 * u1grids(j,i,3) + z1grid(j,i)
        x2h(j,i) = dt/2 * u2grids(j,i,1) + x2grid(j,i)
        y2h(j,i) = dt/2 * u2grids(j,i,2) + y2grid(j,i)
        z2h(j,i) = dt/2 * u2grids(j,i,3) + z2grid(j,i)
     enddo
  enddo

  !!! no need to repara in the half step

  !! obtain geometry property 
  !! drop 1 
  call sphtrans_fwd_real_rsv(nterms,x1mpole,nphi,ntheta,x1h, &
       ctheta,whts,ynms,dwsave)
  call sphtrans_fwd_real_rsv(nterms,y1mpole,nphi,ntheta,y1h, &
       ctheta,whts,ynms,dwsave)    
  call sphtrans_fwd_real_rsv(nterms,z1mpole,nphi,ntheta,z1h, &
       ctheta,whts,ynms,dwsave)    
  !! drop 2
  call sphtrans_fwd_real_rsv(nterms,x2mpole,nphi,ntheta,x2h, &
       ctheta,whts,ynms,dwsave)
  call sphtrans_fwd_real_rsv(nterms,y2mpole,nphi,ntheta,y2h, &
       ctheta,whts,ynms,dwsave)    
  call sphtrans_fwd_real_rsv(nterms,z2mpole,nphi,ntheta,z2h, &
       ctheta,whts,ynms,dwsave)    
 
  call eval_curv_mpole_upAll(1, istep, x1mpole, y1mpole, z1mpole,  &
   h1grid, ws1grid, n1grid, h1mpole, ws1mpole, n1mpole)

  call eval_curv_mpole_upAll(2, istep, x2mpole, y2mpole, z2mpole,  &
   h2grid, ws2grid, n2grid, h2mpole, ws2mpole, n2mpole)

  !! evaluate half step velocity 

   call eval_velo_pair_cond_wmod(istep,   &
   x1h, y1h, z1h, x1mpole, y1mpole, z1mpole,   &
   h1mpole, ws1mpole, n1mpole,    &
   h1grid, ws1grid, n1grid, u1h, en1mpole,    & 
   x2h, y2h, z2h, x2mpole, y2mpole, z2mpole,   &
   h2mpole, ws2mpole, n2mpole,   &
   h2grid, ws2grid, n2grid, u2h, en2mpole)
  
  write(302,*) "istep=",istep,"run=",num_run, "dt=",dt
  do i = 1, ntheta
    do j = 1, nphi
   write(302,*) i,j, u1h(j,i,1), u1h(j,i,2), u1h(j,i,3)
    enddo
  enddo

  !! march using obtained half step velocity, midpoint rule
  do j = 1,nphi
     do i = 1,ntheta
         x1m(j,i) = dt * u1h(j,i,1) + x1grid(j,i)
         y1m(j,i) = dt * u1h(j,i,2) + y1grid(j,i)
         z1m(j,i) = dt * u1h(j,i,3) + z1grid(j,i)
         x2m(j,i) = dt * u2h(j,i,1) + x2grid(j,i)
         y2m(j,i) = dt * u2h(j,i,2) + y2grid(j,i)
         z2m(j,i) = dt * u2h(j,i,3) + z2grid(j,i)
     enddo
  enddo

!!!!!!!! adaptive step !!!!!!!
  max1 = -1.0d0
  max2 = -1.0d0
  call estm_error(ntheta,nphi,march_tol,dt,x1m,y1m,z1m,x1e,y1e,z1e,max1,max_i1,max_j1,newDt1)
  write(301,*) "istep=",istep,"drop1 max rel error=",max1,"at i,j",max_i1,max_j1,"newDt=",newDt1
  call estm_error(ntheta,nphi,march_tol,dt,x2m,y2m,z2m,x2e,y2e,z2e,max2,max_i2,max_j2,newDt2)
  write(301,*) "istep=",istep,"drop2 max rel error=",max2,"at i,j",max_i2,max_j2,"newDt=",newDt2

  if (max1 <= march_tol .and. max2 <= march_tol) then
    dtShrinked = 0  !! flag
  endif 

  !! update dt  
  if (newDt1 < newDt2) then
    dt = newDt1
  else
    dt = newDt2
  endif

  write(301,*) "new dt = ", dt
  num_run = num_run + 1

enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! final update
  do j = 1,nphi
    do i = 1,ntheta
        x1grid(j,i) = x1m(j,i)
        y1grid(j,i) = y1m(j,i)
        z1grid(j,i) = z1m(j,i)
        x2grid(j,i) = x2m(j,i)
        y2grid(j,i) = y2m(j,i)
        z2grid(j,i) = z2m(j,i)
    enddo
 enddo
END SUBROUTINE


SUBROUTINE estm_error(ntheta,nphi,tol,dt,xm,ym,zm,xe,ye,ze,max,max_i,max_j,newDt)

integer :: ntheta, nphi, max_j, max_i
real *8 :: tol, dt, max, newDt
real *8 :: xm(nphi,ntheta), ym(nphi,ntheta), zm(nphi,ntheta), &
  xe(nphi,ntheta), ye(nphi,ntheta), ze(nphi,ntheta)
!! local
real *8,ALLOCATABLE :: error(:,:), ref(:,:)

allocate(error(nphi,ntheta), ref(nphi,ntheta))

do j = 1, nphi
  do i = 1, ntheta
    error(j,i) = (xm(j,i) - xe(j,i)) * (xm(j,i) - xe(j,i))    &
      + (ym(j,i) - ye(j,i)) * (ym(j,i) - ye(j,i))             &
      + (zm(j,i) - ze(j,i)) * (zm(j,i) - ze(j,i)) 
    ref(j,i) =  xm(j,i)*xm(j,i) + ym(j,i)*ym(j,i) + zm(j,i)*zm(j,i)
    error(j,i) = dsqrt(error(j,i) / ref(j,i))
    if (max < error(j,i)) then
      max = error(j,i)
      max_j = j
      max_i = i
    endif
  enddo
enddo

newDt = dt * dsqrt(0.9*tol/max) !! dt can be incremented 

END subroutine
