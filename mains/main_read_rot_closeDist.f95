program main
      
    use mod_assign      
    use mod_upsamp
    use mod_eqsamp
    use mod_orgHelper
    use mod_rotsphf
    use OMP_LIB
    
    implicit real *8 (a-h,o-z)
    
!!      integer ::  nphi, ntheta
!!      integer ::  nterms,mterms, mtheta, mphi
    integer :: iterms, itheta, iphi
!! drop 1
    real *8, ALLOCATABLE :: x1grid(:,:),y1grid(:,:),z1grid(:,:)
    real *8, ALLOCATABLE :: u1grids(:,:,:)      
    complex *16, ALLOCATABLE:: x1mpole(:,:), y1mpole(:,:), z1mpole(:,:)     
    real *8, ALLOCATABLE:: h1grid(:,:), w1grid(:,:), ws1grid(:,:), n1grid(:,:,:)
    complex *16, ALLOCATABLE:: h1mpole(:,:),w1mpole(:,:),ws1mpole(:,:),n1mpole(:,:,:)
    complex *16, ALLOCATABLE:: en1mpole(:,:)
    real *8 :: cen_mass1(3)

!! drop 2
    real *8, ALLOCATABLE :: x2grid(:,:),y2grid(:,:),z2grid(:,:)
    real *8, ALLOCATABLE :: u2grids(:,:,:)    
    complex *16, ALLOCATABLE:: x2mpole(:,:), y2mpole(:,:), z2mpole(:,:)     
    real *8, ALLOCATABLE:: h2grid(:,:), w2grid(:,:), ws2grid(:,:), n2grid(:,:,:)
    complex *16, ALLOCATABLE:: h2mpole(:,:), w2mpole(:,:),ws2mpole(:,:),n2mpole(:,:,:)
    complex *16, ALLOCATABLE:: en2mpole(:,:)
    real *8, ALLOCATABLE:: en2gridsamp(:,:)
    real *8 :: cen_mass2(3)
    
    real *8 :: START,END
    complex *16, ALLOCATABLE:: sphf_s(:,:)
!!! rotate grid
    real *8, ALLOCATABLE:: beta(:),rotmat_one(:,:,:,:)
    real *8, ALLOCATABLE:: xr1grids(:,:,:,:),yr1grids(:,:,:,:),zr1grids(:,:,:,:)
    real *8, ALLOCATABLE:: xr2grids(:,:,:,:),yr2grids(:,:,:,:),zr2grids(:,:,:,:)
    complex *16, ALLOCATABLE:: rx1mpoles(:,:,:,:),ry1mpoles(:,:,:,:),rz1mpoles(:,:,:,:)
    complex *16, ALLOCATABLE:: rx2mpoles(:,:,:,:),ry2mpoles(:,:,:,:),rz2mpoles(:,:,:,:)
! 

    PI = 4.D0*DATAN(1.D0)
    
    center_dist = 2.5d0
    nstep = 2
    
    tmax = 1d0 + 1d-5
    dt = 5d-3  !! initial dt
    
    Ca = 0.1d0
    alpha = 0 ! pi/6
    hintpol = pi/(nterms + 1)/2
    nintpol = 8 + 1
    h_far = 4 * pi/(nterms + 1)

    march_tol = 1d-8

    ratio_cutoff = 1d-2
    dtau = 1d0
    i_repara_max = 20 
    
    call prog_log(6,dt) 
    call prog_log(101,dt)

    call mod_upsamp_mterms_init()
    call grid_dim_init(nterms, ntheta, nphi)    
    call mod_eqsamp_nsamp_init()
    call mod_orgHelper_init()
    call mod_rotsphf_init()
    call mod_upsamp_org_grid_init(theta)
!      
    !! drop 1
    allocate(x1grid(nphi,ntheta),y1grid(nphi,ntheta),z1grid(nphi,ntheta))
    allocate(u1grids(nphi,ntheta,3))
    allocate(h1grid(nphi,ntheta))
    allocate(w1grid(nphi,ntheta))
    allocate(n1grid(nphi,ntheta,3))
    allocate(x1mpole(0:nterms,0:nterms))
    allocate(y1mpole(0:nterms,0:nterms))
    allocate(z1mpole(0:nterms,0:nterms))
    allocate(h1mpole(0:nterms,0:nterms))
    allocate(w1mpole(0:nterms,0:nterms))
    allocate(n1mpole(0:nterms,0:nterms,3))
    allocate(ws1grid(nphi,ntheta))
    allocate(ws1mpole(0:nterms,0:nterms))
    !! full expansion coeff for Enormal 
    allocate(en1mpole(0:nterms,-nterms:nterms))

    !! drop 2 
    allocate(x2grid(nphi,ntheta),y2grid(nphi,ntheta),z2grid(nphi,ntheta))
    allocate(u2grids(nphi,ntheta,3))
    allocate(h2grid(nphi,ntheta))
    allocate(w2grid(nphi,ntheta))
    allocate(n2grid(nphi,ntheta,3))
    allocate(x2mpole(0:nterms,0:nterms))
    allocate(y2mpole(0:nterms,0:nterms))
    allocate(z2mpole(0:nterms,0:nterms))
    allocate(h2mpole(0:nterms,0:nterms))
    allocate(w2mpole(0:nterms,0:nterms))
    allocate(n2mpole(0:nterms,0:nterms,3))
    allocate(ws2grid(nphi,ntheta))
    allocate(ws2mpole(0:nterms,0:nterms))
    !! full expansion coeff for Enormal 
    allocate(en2mpole(0:nterms,-nterms:nterms))

!!ccccccccccccccccccccccccccccccccccccccccccccccc

    !!! read initial grid from file
    
    open(15, file = 'fortcopy151.txt', status='old')
    read(15,*)   !! skip first line
    write(*,*) "read intial grid drop1"
    do i=1,ntheta
      do j=1,nphi
        read(15,*) na, nb, x1grid(j,i), y1grid(j,i), z1grid(j,i) 
        write(*,*) na, nb, x1grid(j,i), y1grid(j,i), z1grid(j,i) 
      enddo
    enddo
    
    close(15)
    
    open(15, file = 'fortcopy251.txt', status='old')
    read(15,*)  !! skip first line
    write(*,*) "read intial grid drop1"
    do i=1,ntheta
      do j=1,nphi
        read(15,*) na, nb, x2grid(j,i), y2grid(j,i), z2grid(j,i) 
        write(*,*) na, nb, x2grid(j,i), y2grid(j,i), z2grid(j,i) 
      enddo
    enddo
    
    close(15)      
    
!!ccccccccccccccccccccccccccccccccccccccccccccccc
!! drop1
    call sphtrans_fwd_real_rsv(nterms,x1mpole,nphi,ntheta,x1grid, &
         ctheta,whts,ynms,dwsave)
    call sphtrans_fwd_real_rsv(nterms,y1mpole,nphi,ntheta,y1grid, &
         ctheta,whts,ynms,dwsave)    
    call sphtrans_fwd_real_rsv(nterms,z1mpole,nphi,ntheta,z1grid, &
         ctheta,whts,ynms,dwsave)    
!! drop 2
    call sphtrans_fwd_real_rsv(nterms,x2mpole,nphi,ntheta,x2grid, &
         ctheta,whts,ynms,dwsave)
    call sphtrans_fwd_real_rsv(nterms,y2mpole,nphi,ntheta,y2grid, &
         ctheta,whts,ynms,dwsave)    
    call sphtrans_fwd_real_rsv(nterms,z2mpole,nphi,ntheta,z2grid, &
         ctheta,whts,ynms,dwsave)    

!!// find closest point to origin on drop 1           
    x0 = 0d0
    y0 = 0d0
    z0 = 0d0  
    allocate(sphf_s(0:nterms,-nterms:nterms))
    call close_dist_flag(0, x0, y0, z0,                  &
    x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole,   &
    dist_cl, theta_cl, phi_cl, sphf_s, xs, ys, zs)
    
    write(*,*) "dist_cl, theta_cl, phi_cl, xs, ys, zs"
    write(*,*) dist_cl, theta_cl, phi_cl, xs, ys, zs

!!//  rotate the grid so that the pole coincides with the closest point to origin

     iphi_cl = nphi/2 + 1
     call rot_grid_pole(theta_cl, iphi_cl, x1mpole, y1mpole, z1mpole,  &
       x1grid, y1grid, z1grid)
     iphi_cl = 1
     call rot_grid_pole(pi-theta_cl,  iphi_cl, x2mpole, y2mpole, z2mpole,  &
       x2grid, y2grid, z2grid)    

    write(199,*) "starting rotated grid: x y z1grid on surface" 
    call log_xyz_grid(199, ntheta, nphi, x1grid, y1grid, z1grid)
    write(299,*) "starting rotated grid: x y z2grid on surface" 
    call log_xyz_grid(299, ntheta, nphi, x2grid, y2grid, z2grid) 

    write(198,*) "starting rotated grid: x y z1mpole on surface" 
    write(298,*) "starting rotated grid: x y z2mpole on surface"       
    do n = 0, nterms
      do m = 0, n
         write(198,*) n,m,x1mpole(n,m), y1mpole(n,m), z1mpole(n,m)
         write(298,*) n,m,x2mpole(n,m), y2mpole(n,m), z2mpole(n,m)
      enddo 
    enddo        

    END program

