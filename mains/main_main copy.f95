      program main
      
      use mod_assign      
      use mod_upsamp
      use mod_eqsamp
      use mod_orgHelper
      use mod_rotsphf
      
      implicit real *8 (a-h,o-z)
      
!!      integer ::  nphi, ntheta
!!      integer ::  nterms,mterms, mtheta, mphi
      integer :: iterms, itheta, iphi

      real *8, ALLOCATABLE:: xgrid(:,:),ygrid(:,:),zgrid(:,:)
      real *8, ALLOCATABLE :: ugrids(:,:,:)
      complex *16, ALLOCATABLE:: xmpole(:,:), ympole(:,:), zmpole(:,:)
      
      complex *16, ALLOCATABLE:: hmpole(:,:),wmpole(:,:),nmpole(:,:,:)
      real *8, ALLOCATABLE:: hgrid(:,:), wgrid(:,:), ngrid(:,:,:)
      
      complex *16, ALLOCATABLE:: mpole(:,:), fgrid(:,:,:,:)
      
      complex *16, ALLOCATABLE:: enmpole(:,:)
      real *8, ALLOCATABLE:: eengrid(:,:)

!ccccccccc    

     call mod_upsamp_mterms_init()
     
!      do n = 0, mterms
!        do m = -n, n
!           do i = 1, mtheta
!              do j = 1, mphi
!                  write(1001,*)  maxup,j,i,n,m,up_sphfgrid(maxup,j,i,n,m),up_d2ydth2(maxup,j,i,n,m)
!              enddo
!            enddo
!          enddo
!        enddo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      call grid_dim_init(nterms, ntheta, nphi)
      
   iup = 1
   iterms = nterms
   itheta = ntheta
   iphi = nphi
   do i = 1, itheta
     do j = 1, iphi
      do n = 0, iterms
        do m = -n, n
                  write(911,*)  j,i,n,m,up_sphfgrid(j,i,n,m,iup), &
            up_dydph(j,i,n,m,iup), up_dydth(j,i,n,m,iup), up_d2ydth2(j,i,n,m,iup), &
            up_d2ydph2(j,i,n,m,iup), up_d2ydthdph(j,i,n,m,iup)
              enddo
            enddo
          enddo
        enddo
      
      call mod_eqsamp_nsamp_init()
      
!      do n = 0, nterms
!        do m = -n, n
!           do i = 1, nsamp
!              do j = 1, nphi
!                  write(1002,*)  j,i,n,m,eyfgrid(j,i,n,m)
!              enddo
!            enddo
!          enddo
!        enddo

!
!      allocate(mpole(0:nterms,-nterms:nterms))
!!      allocate(fgrid(nphi,ntheta,0:nterms,-nterms:nterms))
!
!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call mod_orgHelper_init()
      call mod_rotsphf_init()
      call mod_upsamp_org_grid_init(theta)
      
!      do n = 0, nterms
!        do m = -n, n
!           do i = 1, ntheta
!              do j = 1, nphi
!                  write(1004,*)  j,i,n,m,sphfgrid(j,i,n,m)
!              enddo
!            enddo
!          enddo
!        enddo
!        
!
!       do n = 0, nterms
!         do m = -n, n
!          do l = 1, ntheta
!          do k = 1, nphi
!           do i = 1, ntheta
!             do j = 1, nphi
!                write(1006,*)  j,i,k,l,n,m,rsphfgrids(j,i,k,l,n,m)
!             enddo
!           enddo
!           enddo
!           enddo
!         enddo
!       enddo
!
!      do n = 0, mterms
!        do m = -n, n
!           do i = 1, ntheta
!              do j = 1, nphi
!                  write(1008,*)  j,i,n,m,up_yfngrid(j,i,n,m)
!              enddo
!            enddo
!          enddo
!        enddo

      allocate(xgrid(nphi,ntheta),ygrid(nphi,ntheta),zgrid(nphi,ntheta))
      allocate(ugrids(nphi,ntheta,3))
      allocate(hgrid(nphi,ntheta))
      allocate(wgrid(nphi,ntheta))
      allocate(ngrid(nphi,ntheta,3))
     
      allocate(xmpole(0:nterms,0:nterms))
      allocate(ympole(0:nterms,0:nterms))
      allocate(zmpole(0:nterms,0:nterms))
      allocate(hmpole(0:nterms,0:nterms))
      allocate(wmpole(0:nterms,0:nterms))
      allocate(nmpole(0:nterms,0:nterms,3))
      !! full expansion coeff for Enormal 
      allocate(enmpole(0:nterms,-nterms:nterms))
      allocate(eengrid(nphi,nsamp))
      
      !!! initial grid on a unit sphere
      !!! build 2 parameter grids phigrid and thetagrid
      PI = 4.D0*DATAN(1.D0)
      phi = 2*PI/nphi
      do i=1,ntheta
        write(*,*) i,dsin(theta(i))
        do j=1,nphi
          phij = (j-1)*phi
          xgrid(j, i) = dcos(phij)*dsin(theta(i)) !* 0.99d0
          ygrid(j, i) = dsin(phij)*dsin(theta(i)) !* 0.99d0
          zgrid(j, i) = ctheta(i) !
        enddo
      enddo
      write(52,*) "starting grid: x y zgrid on sphere" 
      call log_xyz_grid(52, ntheta, nphi, xgrid,ygrid,zgrid)
      
!!ccccccccccccccccccccccccccccccccccccccccccccccc

      call sphtrans_fwd_real_rsv(nterms,xmpole,nphi,ntheta,xgrid, &
           ctheta,whts,ynms,dwsave)
      call sphtrans_fwd_real_rsv(nterms,ympole,nphi,ntheta,ygrid, &
           ctheta,whts,ynms,dwsave)    
      call sphtrans_fwd_real_rsv(nterms,zmpole,nphi,ntheta,zgrid, &
           ctheta,whts,ynms,dwsave)    


!
      nstep = 1
      dt = 5d-3
      Ca = 0.1d0
!      
do istep = 1, nstep
!!!ccccccccccc    
      write(14,*) "istep =", istep  
      write(15,*) "istep =", istep
      write(16,*) "istep =", istep 
      write(19,*) "istep =", istep 
      write(20,*) "istep =", istep 
      write(21,*) "istep =", istep 
      write(32,*) "istep =", istep 
      write(88,*) "istep =", istep 
      write(55,*) "istep =", istep 
      write(*,*) "istep =", istep ,"dt =", dt
!!!ccccccccccc      
!!ccccccccccccccccccccccccccccccccccccccccccccccc
!!ccccccc evaluate curvature for x,y,zmpole ccccc   

      call eval_curv_mpole_uphonly(xmpole, ympole, zmpole, &
       hgrid, wgrid, ngrid, hmpole)
      write(*,*) "eval_curv done"
      call sphtrans_real(nterms,hmpole,nphi,ntheta,hgrid, &
        ctheta,ynms,dwsave)   
      write(*,*) "obtain hgrid"      
      do i = 1,ntheta
        do j = 1, nphi
          write(183,*) i,j,"H",hgrid(j,i),"W",wgrid(j,i) 
        enddo
      enddo
      
!!!!!!!!!
      x0 = 1.5d0
      y0 = 0d0
      z0 = -0.8d0

!      call close_dist_test(x0, y0, z0,                   &
!       xgrid, ygrid, zgrid, xmpole, ympole, zmpole, &
!       theta_cl, phi_cl, dist_cl)
!      
!     write(*,*) "closest point:", theta_cl, phi_cl, "dist =", dist_cl
!      !!call nearSing_integral_test(xgrid,ygrid,zgrid,wgrid,ngrid)
!      call nearSing_integral_test_upsamp(xmpole, ympole, zmpole, &
!          xgrid, ygrid, zgrid, hgrid, wgrid, ngrid)      
!     
      write(*,*) "test done"
      call eval_area_vol_real(ntheta, nphi, theta, whts, &
       xgrid, ygrid, zgrid, wgrid, ngrid, surf_area, tot_vol)
      write(*,*) "eval area and vol"
!!ccccccccccccccccccccccccccccccccccccccccccccccc
     
!      call eval_velo(istep, Ca, &
!      xgrid, ygrid, zgrid, xmpole, ympole, zmpole, &
!      hgrid, wgrid, ngrid, &
!      ugrids, enmpole)
       call eval_velo_single_wmod(istep, Ca,  &
       xgrid, ygrid, zgrid, xmpole, ympole, zmpole,  &
       hgrid, wgrid, ngrid, ugrids, enmpole)
       
       write(*,*) "flow obtained"
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      call expl_fwd_elr(ntheta, nphi, dt, ugrids, xgrid, ygrid, zgrid) 
!!cccccccccccc

      nup = 2
      ratio_cutoff = 4d-1
      dtau = 1d0
      i_repara_max = 20 
      
      call repara_real_upsamp(nup, ratio_cutoff, dtau, i_repara_max,  &
       xgrid, ygrid, zgrid, xmpole, ympole, zmpole)
     
      !! transform x,y,zmpole to xgrid, ygrid, zgrid
      call sphtrans_real(nterms,xmpole,nphi,ntheta,  & 
        xgrid,ctheta,ynms,dwsave)
      write(*,*) "*nterms =", nterms, "ntheta =", ntheta, "nphi =",nphi
      call sphtrans_real(nterms,ympole,nphi,ntheta,  & 
        ygrid,ctheta,ynms,dwsave)    
      call sphtrans_real(nterms,zmpole,nphi,ntheta,  & 
        zgrid,ctheta,ynms,dwsave)         
     
      write(52,*) "istep =", istep
      write(52,*) "trunc nextGrid final x y z"  
      call log_xyz_grid(52, ntheta, nphi, xgrid, ygrid, zgrid)

!!      write(53,*) "istep =", istep
!!      write(53,*) "equal angle x y zgrid"  
!!      call log_xyz_grid(53, nsamp, nphi, exgrid, eygrid, ezgrid)
!!      mid = (nsamp+1)/2
!!      shortaxis = exgrid(1,mid)**2 + eygrid(1,mid)**2 + ezgrid(1,mid)**2
!!      shortaxis = dsqrt(shortaxis)
!!      longaxis = ezgrid(1,1)
!!      write(56,*) "istep =",istep,longaxis,shortaxis,longaxis/shortaxis

       call mod_eqsamp_eval_grid_xyz(istep, xmpole, ympole, zmpole)
enddo

      write(101,*) "All prameters"
      write(101,*)  "nup =", nup, "ratio_cutoff =", ratio_cutoff,  & 
        "dtau =", dtau, "i_repara_max =", i_repara_max
      write(101,*)  "nstep =", nstep,   & 
        "dt =", dt, "Ca =", Ca
      write(101,*) "nterms =", nterms, "ntheta =", ntheta,  & 
        "maxup =", maxup 

      END program
