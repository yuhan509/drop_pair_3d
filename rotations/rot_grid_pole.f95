!!! input: 
!!!   x,y,zmpole; 
!!!   new pole position for drop 1 in surface coordinate: theta_cl, iphi_cl; 
!!! output: 
!!!   rotated x,y,zmpole
!!!   rotated x,y,zgrid

SUBROUTINE rot_grid_pole(theta_cl, iphi_cl, xmpole, ympole, zmpole,  &
  xgrid, ygrid, zgrid)

  use mod_orgHelper

  implicit real *8 (a-h,o-z)
  
  real*8,INTENT(IN)  :: theta_cl
  integer,INTENT(IN) :: iphi_cl
  complex*16 :: xmpole(0:nterms,0:nterms), ympole(0:nterms,0:nterms), zmpole(0:nterms,0:nterms)
  real*8 :: xgrid(nphi, ntheta), ygrid(nphi, ntheta), zgrid(nphi, ntheta)

  real *8, ALLOCATABLE:: beta(:), rotmat_one(:,:,:,:)
  complex *16, ALLOCATABLE:: rxmpoles(:,:,:,:),rympoles(:,:,:,:),rzmpoles(:,:,:,:)
  real *8, ALLOCATABLE:: xrgrids(:,:,:,:),yrgrids(:,:,:,:),zrgrids(:,:,:,:)

    nrot = nphi
    nbeta = 1
    allocate(rotmat_one(0:nterms,-nterms:nterms,-nterms:nterms,nbeta))
    allocate(beta(nbeta))
    allocate(xrgrids(nphi,ntheta,nrot,nbeta),yrgrids(nphi,ntheta,nrot,nbeta), &
    zrgrids(nphi,ntheta,nrot,nbeta))
    allocate(rxmpoles(0:nterms,0:nterms,nrot,nbeta), rympoles(0:nterms,0:nterms,nrot,nbeta), &
    rzmpoles(0:nterms,0:nterms,nrot,nbeta))

    beta(1) = theta_cl
    call rotgrid_fsr_real_init(nterms,nbeta,beta,rotmat_one)
    
    call rotgrid_fsr_real_mpolesout(nterms,xmpole,nphi,ntheta,  &
      nbeta,beta,nrot,xrgrids,rotmat_one,ctheta,ynms,dwsave,rxmpoles)
      
    call rotgrid_fsr_real_mpolesout(nterms,ympole,nphi,ntheta,  &
      nbeta,beta,nrot,yrgrids,rotmat_one,ctheta,ynms,dwsave,rympoles)
      
    call rotgrid_fsr_real_mpolesout(nterms,zmpole,nphi,ntheta,  &
      nbeta,beta,nrot,zrgrids,rotmat_one,ctheta,ynms,dwsave,rzmpoles)
      
    k = iphi_cl
    do i = 1, ntheta
      do j = 1, nphi
        xgrid(j,i) = xrgrids(j,i,k,1)
        ygrid(j,i) = yrgrids(j,i,k,1)
        zgrid(j,i) = zrgrids(j,i,k,1)
      enddo
    enddo

    do n = 0, nterms
      do m = 0, n
        xmpole(n,m) = rxmpoles(n,m,k,1)
        ympole(n,m) = rympoles(n,m,k,1)
        zmpole(n,m) = rzmpoles(n,m,k,1)
      enddo 
    enddo  
    
END SUBROUTINE


!!!!// find closest point to origin on drop 1           
!!    x0 = 0d0
!!    y0 = 0d0
!!    z0 = 0d0  
!!    allocate(sphf_s(0:nterms,-nterms:nterms))
!!    call close_dist(x0, y0, z0,                          &
!!    x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole,   &
!!    dist_cl, theta_cl, phi_cl, sphf_s, xs, ys, zs)
!!    
!!    write(*,*) "dist_cl, theta_cl, phi_cl, xs, ys, zs"
!!    write(*,*) dist_cl, theta_cl, phi_cl, xs, ys, zs
!!
!!!!//  rotate the grid so that the pole coincides with the closest point to origin
!!
!!     iphi_cl = nphi/2 + 1
!!     call rot_grid_pole(theta_cl, iphi_cl, x1mpole, y1mpole, z1mpole,  &
!!       x1grid, y1grid, z1grid)
!!     iphi_cl = 1
!!     call rot_grid_pole(pi-theta_cl,  iphi_cl, x2mpole, y2mpole, z2mpole,  &
!!       x2grid, y2grid, z2grid)    
