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
      complex *16, ALLOCATABLE:: x1mpole(:,:), y1mpole(:,:), z1mpole(:,:), sphf_s(:,:)     
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
! 

      PI = 4.D0*DATAN(1.D0)
      
      center_dist = 3d0
      nstep = 100000
      tmax = 100d0 + 1d-5
      
      Ca = 0.1d0
      alpha = pi/4
      hintpol = pi/(nterms + 1)/3
      nintpol = 8 + 1
      h_far = 4 * pi/(nterms + 1)
      dt = 5d-3
      march_tol = 5d-7

      ratio_cutoff = 1d-2
      dtau = 1d0
      i_repara_max = 8 
      call prog_log(101,dt)
      call prog_log(6,dt)
  
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

      allocate(sphf_s(0:nterms, -nterms:nterms))   !! used in close_dist
!!ccccccccccccccccccccccccccccccccccccccccccccccc

!!! 

      !!! initial grid on a unit sphere

      phi = 2*PI/nphi
      c = 1d0 !! (x*x + y*y)/(b*b) + z*z/(c*c) = 1 ; volume = 4/3*pi*b*b*c
      b = 1d0 !dsqrt(1d0/c)
      a = 1d0 
      do i=1,ntheta
        write(*,*) i,dsin(theta(i))
        do j=1,nphi
          phij = (j-1)*phi
          !! drop 1
          x1grid(j, i) = dcos(phij)*dsin(theta(i))*a 
          y1grid(j, i) = dsin(phij)*dsin(theta(i))*b 
          z1grid(j, i) = ctheta(i)*c + center_dist/2
          !! drop 2
          x2grid(j, i) = dcos(phij)*dsin(theta(i))*a
          y2grid(j, i) = dsin(phij)*dsin(theta(i))*b
          z2grid(j, i) = ctheta(i)*c  - center_dist/2
        enddo
      enddo
      write(152,*) "starting grid: x y z1grid on sphere" 
      call log_xyz_grid(152, ntheta, nphi, x1grid, y1grid, z1grid)
      write(252,*) "starting grid: x y z1grid on sphere" 
      call log_xyz_grid(252, ntheta, nphi, x2grid, y2grid, z2grid)   
      
      !!! mean curvature initial ellipsoid
      write(189,*) "mean curvature initial ellipsoid"
      do i =1, ntheta
        tmp = ctheta(i)**2 + c**2 * dsin(theta(i))**2
        write(189,*) i, - 0.5d0 * c*(1d0 + tmp) / dsqrt(tmp)**3
      enddo
      write(189,*) "infinitesimal area element / sin(theta(i))"
      do i =1, ntheta
        write(189,*) i, dsqrt(dsin(theta(i))**2*c**2 + ctheta(i)**2)
      enddo
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

      !! equally distributed theta sample grid for En2
      allocate(en2gridsamp(nphi,nsamp))
!      
   
      time_sum = 0d0
      istep = 0
do while(time_sum < tmax .and. istep < nstep)
!!!ccccccccccc    
      istep = istep + 1
      write(*,*) "**** istep =", istep, "t =",time_sum
      write(300,*) "istep =", istep, "t =", time_sum
      write(50,*) "istep =", istep 
!!ccccccccccccccccccccccccccccccccccccccccccccccc
!!ccccccc evaluate curvature ccccccc   
      !! drop 1
  !    call eval_curv_mpole_uphonly(x1mpole, y1mpole, z1mpole, &
  !     h1grid, w1grid, n1grid, h1mpole)
      call eval_curv_mpole_upAll(1, istep, x1mpole, y1mpole, z1mpole,  &
       h1grid, ws1grid, n1grid, h1mpole, ws1mpole, n1mpole)
      write(*,*) "eval_curv done"
      call sphtrans_real(nterms,h1mpole,nphi,ntheta,h1grid,  &
        ctheta,ynms,dwsave)   
      write(*,*) "obtain h1grid"      
      call eval_area_vol_cenMass_real(ntheta, nphi, whts,  &
      x1grid, y1grid, z1grid, ws1grid, n1grid, surf_area1, tot_vol1, cen_mass1)

      !! drop 2
 !     call eval_curv_mpole_uphonly(x2mpole, y2mpole, z2mpole, &
 !      h2grid, w2grid, n2grid, h2mpole)
       call eval_curv_mpole_upAll(2, istep, x2mpole, y2mpole, z2mpole,  &
       h2grid, ws2grid, n2grid, h2mpole, ws2mpole, n2mpole)
      write(*,*) "eval_curv done"
      call sphtrans_real(nterms,h2mpole,nphi,ntheta,h2grid, &
        ctheta,ynms,dwsave)   
      write(*,*) "obtain h2grid"      
       call eval_area_vol_cenMass_real(ntheta, nphi, whts, &
       x2grid, y2grid, z2grid, ws2grid, n2grid, surf_area2, tot_vol2, cen_mass2)
      
      write(141,*) "istep=",istep-1, surf_area1, tot_vol1
      write(142,*) "istep=",istep-1, cen_mass1(1), cen_mass1(2), cen_mass1(3)
      write(241,*) "istep=",istep-1, surf_area2, tot_vol2
      write(242,*) "istep=",istep-1, cen_mass2(1), cen_mass2(2), cen_mass2(3)
      write(183,*) "istep=",istep-1
      write(283,*) "istep=",istep-1
      do i = 1,ntheta
        do j = 1, nphi
          write(183,*) i,j,"H",h1grid(j,i),"WS",ws1grid(j,i), &
          "nxyz",n1grid(j,i,1), n1grid(j,i,2), n1grid(j,i,3)
          write(283,*) i,j,"H",h2grid(j,i),"WS",ws2grid(j,i), &
          "nxyz",n2grid(j,i,1), n2grid(j,i,2), n2grid(j,i,3)
        enddo
      enddo    
!!ccccccccccccccccccccccccccccccccccccccccccccccc

      call eval_velo_pair_cond_wmod(istep,                    &
      x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole, &
      h1mpole, ws1mpole, n1mpole, &
      h1grid, ws1grid, n1grid, u1grids, en1mpole, & 
      x2grid, y2grid, z2grid, x2mpole, y2mpole, z2mpole, &
      h2mpole, ws2mpole, n2mpole, &
      h2grid, ws2grid, n2grid, u2grids, en2mpole)
      
      
          write(134,*) "istep=",istep,"u1grid x,y,z component"
          write(234,*) "istep=",istep,"u2grid x,y,z component"
          do i = 1, ntheta
            do j = 1, nphi
              write(134,*) i,j, u1grids(j,i,1), u1grids(j,i,2), u1grids(j,i,3)
              write(234,*) i,j, u2grids(j,i,1), u2grids(j,i,2), u2grids(j,i,3)
            enddo
          enddo
      

!      call mod_eqsamp_eval_grid_fcomplex(en2mpole, en2gridsamp)
       call mod_eqsamp_eval_grid_freal(en2mpole, en2gridsamp)
      do i = 1,nsamp
       do j = 1, nphi
          write(110,*) j,i,en2gridsamp(j,i)
       enddo
      enddo
         
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!      call expl_fwd_elr(ntheta, nphi, dt, u1grids, x1grid, y1grid, z1grid)
!       
!      call expl_fwd_elr(ntheta, nphi, dt, u2grids, x2grid, y2grid, z2grid) 

      !! adaptive time step could change dt
      call  expl_midp_pair(istep, dt,         &
        u1grids, x1grid, y1grid, z1grid,      &
        u2grids, x2grid, y2grid, z2grid) 
!!cccccccccccc

      write(1510,*) "istep =", istep, " after marching before repara  x y z, drop1"  
      call log_xyz_grid(1510, ntheta, nphi, x1grid, y1grid, z1grid)
      
      START = omp_get_wtime()

      call repara_real_upsamp(1, istep, x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole)
       END = omp_get_wtime() 

      !! transform x,y,z1mpole to x1grid, y1grid, z1grid
      call sphtrans_real(nterms,x1mpole,nphi,ntheta,  & 
        x1grid,ctheta,ynms,dwsave)
      call sphtrans_real(nterms,y1mpole,nphi,ntheta,  & 
        y1grid,ctheta,ynms,dwsave)    
      call sphtrans_real(nterms,z1mpole,nphi,ntheta,  & 
        z1grid,ctheta,ynms,dwsave)     

      write(151,*) "istep =", istep, "after repara  x y z, drop1"  
      call log_xyz_grid(151, ntheta, nphi, x1grid, y1grid, z1grid)
      call log_xyz_grid(152, ntheta, nphi, x1grid, y1grid, z1grid)
      
      call repara_real_upsamp(2, istep, x2grid, y2grid, z2grid, x2mpole, y2mpole, z2mpole)
     
      call sphtrans_real(nterms,x2mpole,nphi,ntheta,  & 
        x2grid,ctheta,ynms,dwsave)
      call sphtrans_real(nterms,y2mpole,nphi,ntheta,  & 
        y2grid,ctheta,ynms,dwsave)    
      call sphtrans_real(nterms,z2mpole,nphi,ntheta,  & 
        z2grid,ctheta,ynms,dwsave)         
     
        PRINT *, "total repara work took", END - START, "seconds"    
        START = omp_get_wtime()

      write(251,*) "istep =", istep, "after repara  x y z, drop2"  
      call log_xyz_grid(251, ntheta, nphi, x2grid, y2grid, z2grid)
      call log_xyz_grid(252, ntheta, nphi, x2grid, y2grid, z2grid)

      call mod_eqsamp_eval_grid_xyz_pair(istep, &
      x1mpole, y1mpole, z1mpole, x2mpole, y2mpole, z2mpole)
      
      write(*,*) "dt =", dt
      time_sum = time_sum + dt
!! axisymmetric case close distance to origin     
!      if (ez1grid(nphi,nsamp) .lt. 5d-2) then
!        write(101,*) "stop at ez1grid(nphi,nsamp) .lt. 5d-2"
!        exit
!      endif
!!ccccccccccccccc

!! non-axisymmetric case close distance to origin
      dist_min = 1d3
      x0 = 0d0
      y0 = 0d0
      z0 = 0d0      
      do i = 1, ntheta
         do j = 1, nphi
            dist = (x0 - x1grid(j,i))**2 + (y0 - y1grid(j,i))**2  + (z0 - z1grid(j,i))**2
            if (dist_min .gt. dist) then
              dist_min = dist
              j_cl = j
              i_cl = i
            endif
         enddo
      enddo
      dist_min = dsqrt(dist_min)
     !! non-axisymmetric case close distance to origin
      write(102,*)  "istep =",istep, "dist_min =",dist_min, "i,j",i_cl, j_cl
      write(*,*)  "istep =",istep, "dist_min =",dist_min, "i,j",i_cl, j_cl
      if (dist_min .lt. 5d-2) then
        write(101,*) "stop at dist_min .lt. 5d-2"
        exit
      endif
!!cccccccccccccc 

enddo



      END program
