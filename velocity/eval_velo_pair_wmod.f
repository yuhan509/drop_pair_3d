      SUBROUTINE eval_velo_pair_wmod(istep, Ca, 
     $ x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole, 
     $ h1grid, w1grid, n1grid, u1grid, en1mpole,  
     $ x2grid, y2grid, z2grid, x2mpole, y2mpole, z2mpole, 
     $ h2grid, w2grid, n2grid, u2grid, en2mpole)
    
      use mod_orgHelper
      use mod_rotsphf
      use mod_upsamp
      
      implicit real *8 (a-h,o-z)
      
      integer :: istat = 2, istep
      integer :: nintpol = 8
      real *8 :: Ca

      real *8 :: x1grid(nphi,ntheta),y1grid(nphi,ntheta),
     $ z1grid(nphi,ntheta),h1grid(nphi,ntheta),w1grid(nphi,ntheta), 
     $ n1grid(nphi,ntheta,3)
      complex *16 :: x1mpole(0:nterms,0:nterms),
     $ y1mpole(0:nterms,0:nterms),z1mpole(0:nterms,0:nterms)

      real *8 :: x2grid(nphi,ntheta),y2grid(nphi,ntheta),
     $ z2grid(nphi,ntheta),h2grid(nphi,ntheta),w2grid(nphi,ntheta), 
     $ n2grid(nphi,ntheta,3)
      complex *16 :: x2mpole(0:nterms,0:nterms),
     $ y2mpole(0:nterms,0:nterms),z2mpole(0:nterms,0:nterms)

      !! output 
      real *8 :: u1grid(nphi,ntheta,3), u2grid(nphi,ntheta,3)
      complex *16 :: en1mpole(0:nterms,-nterms:nterms),
     $ en2mpole(0:nterms,-nterms:nterms)
      !! end of output
      
      complex *16, ALLOCATABLE :: sphf_s(:,:)
      integer :: iterms, itheta, iphi
      real *8, ALLOCATABLE :: whts_i(:)
      real *8, dimension(:,:),ALLOCATABLE :: 
     $ x1igrid, y1igrid, z1igrid, ws1igrid, h1igrid,
     $ x2igrid, y2igrid, z2igrid, ws2igrid, h2igrid
      real *8, dimension(:,:,:),ALLOCATABLE :: n1igrid, n2igrid
      complex *16, ALLOCATABLE :: sphfigrid(:,:,:,:)

      complex *16,dimension(:,:),ALLOCATABLE:: cEinfn1grid, cEinfn2grid
      real *8, dimension(:,:), ALLOCATABLE:: ws1grid, ws2grid
      complex *16, dimension(:,:), ALLOCATABLE:: ws1mpole, ws2mpole
      real *8, dimension(:,:,:,:), ALLOCATABLE:: rx1grids, ry1grids,
     $ rz1grids, rws1grids, rx2grids, ry2grids, rz2grids, rws2grids

      complex *16,dimension(:,:,:,:),ALLOCATABLE :: gxs1d1grids,
     $ gys1d1grids, gzs1d1grids, gxs2d2grids, gys2d2grids, gzs2d2grids
      complex *16,dimension(:,:,:,:),ALLOCATABLE :: gx1mpoles,gx2mpoles,
     $ gy1mpoles,gy2mpoles,gz1mpoles,gz2mpoles

      complex *16,dimension(:,:,:,:),ALLOCATABLE :: f1mpoles,f2mpoles

      complex *16, dimension(:,:,:,:), ALLOCATABLE :: cls1d1grids,
     $ cls1d2grids, cls2d1grids, cls2d2grids
      complex *16, dimension(:,:,:,:), ALLOCATABLE :: cls1d1mpoles,
     $ cls1d2mpoles, cls2d1mpoles, cls2d2mpoles

      complex *16, dimension(:,:), ALLOCATABLE:: Einfn1mpole,Einfn2mpole

      complex *16, ALLOCATABLE:: cefield_mat(:,:), cEinf_col(:) 

      INTEGER :: NRHS=1,N2
      INTEGER :: INFO
      INTEGER, allocatable :: ipiv(:)

      real *8, ALLOCATABLE, dimension(:,:):: En1, En2
      real *8, dimension(:,:), allocatable:: fx1,fy1,fz1,fx2,fy2,fz2
      complex *16, dimension(:,:), allocatable :: fx1mpole, fy1mpole,
     $ fz1mpole, fx2mpole, fy2mpole, fz2mpole
      real *8, dimension(:,:,:,:), allocatable :: rfx1grids, rfy1grids,
     $ rfz1grids, rfx2grids, rfy2grids, rfz2grids

      real *8, dimension(:,:,:), allocatable :: us1d1grid, us2d1grid, 
     $ us1d2grid, us2d2grid
      complex *16, dimension(:,:,:), allocatable :: u1mpoles, u2mpoles
      real *8, dimension(:,:), allocatable :: fx1igrid, fy1igrid,
     $ fz1igrid, fx2igrid, fy2igrid, fz2igrid, En1igrid, En2igrid
      
      complex* 16::ctmp,cquadr,ctmp1,ctmp2,ctmp3,cquadr1,cquadr2,cquadr3
     $ , cres
      !! dimension of matrix block (exclude n=0,m=0 in expansions for En1, En2)
      N2=((nterms+1)*(nterms+1) - 1)*2

      PI = 4.D0*DATAN(1.D0)
cccccccccccccccccccccccccccccccccccccccccccccc
      allocate(ws1grid(nphi,ntheta), ws2grid(nphi,ntheta))
      allocate(cEinfn1grid(nphi,ntheta),cEinfn2grid(nphi,ntheta))
      
      alpha = pi/6
      
      do i=1,ntheta
        do j=1,nphi
          ws1grid(j,i) = w1grid(j,i)/dsin(theta(i))
          ws2grid(j,i) = w2grid(j,i)/dsin(theta(i))
          !! E_inf * ngrid
          cEinfn1grid(j,i) = dcmplx(n1grid(j,i,3)*dcos(alpha) + 
     $     n1grid(j,i,1)*dsin(alpha), 0d0)
          cEinfn2grid(j,i) = dcmplx(n2grid(j,i,3)*dcos(alpha) + 
     $     n2grid(j,i,1)*dsin(alpha), 0d0)
          !write(88,*) j, i, ngrid(j,i,3), ctheta(i)
        enddo
      enddo
      write(*,*) "Einfn assign done"

      allocate(ws1mpole(0:nterms,0:nterms),ws2mpole(0:nterms,0:nterms))

      call sphtrans_fwd_real_rsv(nterms,ws1mpole,nphi,ntheta,ws1grid,
     $  ctheta,whts,ynms,dwsave)
      call sphtrans_fwd_real_rsv(nterms,ws2mpole,nphi,ntheta,ws2grid,
     $  ctheta,whts,ynms,dwsave)
cccccccccccccccccccccccccccccccccccc
cccccc rotate x,y,z mpole and local area element ccccccccccc
      nbeta = ntheta
      nrot = nphi

      allocate(rx1grids(nphi,ntheta,nrot,nbeta),
     $         ry1grids(nphi,ntheta,nrot,nbeta),
     $         rz1grids(nphi,ntheta,nrot,nbeta),
     $         rws1grids(nphi,ntheta,nrot,nbeta))

      allocate(rx2grids(nphi,ntheta,nrot,nbeta),
     $         ry2grids(nphi,ntheta,nrot,nbeta),
     $         rz2grids(nphi,ntheta,nrot,nbeta),
     $         rws2grids(nphi,ntheta,nrot,nbeta))

      call rotgrid_dsr_real(nterms,x1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rx1grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,y1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,ry1grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,z1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rz1grids,
     $ rotmat,ctheta,ynms,dwsave)     
      call rotgrid_dsr_real(nterms,ws1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rws1grids,
     $ rotmat,ctheta,ynms,dwsave)     

      call rotgrid_dsr_real(nterms,x2mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rx2grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,y2mpole,nphi,ntheta,
     $ nbeta,theta,nrot,ry2grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,z2mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rz2grids,
     $ rotmat,ctheta,ynms,dwsave)     
      call rotgrid_dsr_real(nterms,ws2mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rws2grids,
     $ rotmat,ctheta,ynms,dwsave)     


      allocate(cls1d1grids(nrot,nbeta,0:nterms,-nterms:nterms),
     $  cls2d2grids(nrot,nbeta,0:nterms,-nterms:nterms), 
     $  cls1d2grids(nrot,nbeta,0:nterms,-nterms:nterms), 
     $  cls2d1grids(nrot,nbeta,0:nterms,-nterms:nterms))

      allocate(gxs1d1grids(nrot,nbeta,0:nterms,-nterms:nterms),
     $ gys1d1grids(nrot,nbeta,0:nterms,-nterms:nterms),
     $ gzs1d1grids(nrot,nbeta,0:nterms,-nterms:nterms),
     $ gxs2d2grids(nrot,nbeta,0:nterms,-nterms:nterms),
     $ gys2d2grids(nrot,nbeta,0:nterms,-nterms:nterms),
     $ gzs2d2grids(nrot,nbeta,0:nterms,-nterms:nterms))

      write(*,*) "start: special quadr efield"
cc    special gauss-legendre quadrature(lapalce) for all rotated grids
      do n = 0,nterms
        do m = -n,n
          !!write(*,*) "using rotating sph", n, m
          do kk = 1,nrot
            do k = 1,nbeta
              write(*,*) "s1d1 n,m,kk,k", n,m,kk,k
              call integral_eval_efield_cmpl3(nterms,ntheta,nphi,
     $        rx1grids(:,:,kk,k),ry1grids(:,:,kk,k),rz1grids(:,:,kk,k),
     $        rws1grids(:,:,kk,k), rsphfgrids(:,:,kk,k,n,m), swhts,
     $        x1grid(kk,k), y1grid(kk,k), z1grid(kk,k), 
     $        cquadr1, cquadr2, cquadr3)  
              !! naming : c- is for complex; l- is for laplace kernel;
              !! s1- means source(kernel) point is on drop 1;
              !! d1- means the contribution from integrating over drop 1.
              !! (kk,k) represents kernel point, will serve as 
              cls1d1grids(kk,k,n,m) = cquadr1 * n1grid(kk,k,1)
     $         + cquadr2 * n1grid(kk,k,2) + cquadr3 * n1grid(kk,k,3)
              gxs1d1grids(kk, k, n, m) = cquadr1
              gys1d1grids(kk, k, n, m) = cquadr2
              gzs1d1grids(kk, k, n, m) = cquadr3
              write(9411,*) kk,k,n,m,cls1d1grids(kk,k,n,m)
            enddo
          enddo
        enddo
      enddo     

cc    special gauss-legendre quadrature(lapalce) for all rotated grids
      
      do n = 0,nterms
        do m = -n,n
          !!write(*,*) "using rotating sph", n, m
          do kk = 1,nrot
            do k = 1,nbeta
              write(*,*) "s2d2 n,m,kk,k", n,m,kk,k
              call integral_eval_efield_cmpl3(nterms,ntheta,nphi,
     $        rx2grids(:,:,kk,k),ry2grids(:,:,kk,k),rz2grids(:,:,kk,k),
     $        rws2grids(:,:,kk,k), rsphfgrids(:,:,kk,k,n,m), swhts,
     $        x2grid(kk,k), y2grid(kk,k), z2grid(kk,k), 
     $        cquadr1, cquadr2, cquadr3)  

              cls2d2grids(kk,k,n,m) = cquadr1 * n2grid(kk,k,1)
     $         + cquadr2 * n2grid(kk,k,2) + cquadr3 * n2grid(kk,k,3)
              gxs2d2grids(kk, k, n, m) = cquadr1
              gys2d2grids(kk, k, n, m) = cquadr2
              gzs2d2grids(kk, k, n, m) = cquadr3
              write(9422,*) kk,k,n,m,cls2d2grids(kk,k,n,m)
            enddo
          enddo
        enddo
      enddo     

cc    prepare expansion of special quadratures for interpolation 
!      allocate(
!     $ f1mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms),
!     $ f2mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms))
!      do n = 0, nterms
!        do m = -n, n
!          call sphtrans_fwd_cmpl_rsv(nterms,f1mpoles(:,:,n,m),
!     $    nphi,ntheta,cls1d1grids(:,:,n,m),ctheta,whts,ynms,zwsave)
!          call sphtrans_fwd_cmpl_rsv(nterms,f2mpoles(:,:,n,m),
!     $    nphi,ntheta,cls2d2grids(:,:,n,m),ctheta,whts,ynms,zwsave)
!        enddo
!      enddo

      allocate(
     $ gx1mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms),
     $ gy1mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms),
     $ gz1mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms),
     $ gx2mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms),
     $ gy2mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms),
     $ gz2mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms))
      do n = 0, nterms
        do m = -n, n
          call sphtrans_fwd_cmpl_rsv(nterms,gx1mpoles(:,:,n,m),
     $    nphi,ntheta,gxs1d1grids(:,:,n,m),ctheta,whts,ynms,zwsave)
          call sphtrans_fwd_cmpl_rsv(nterms,gy1mpoles(:,:,n,m),
     $    nphi,ntheta,gys1d1grids(:,:,n,m),ctheta,whts,ynms,zwsave)
          call sphtrans_fwd_cmpl_rsv(nterms,gz1mpoles(:,:,n,m),
     $    nphi,ntheta,gzs1d1grids(:,:,n,m),ctheta,whts,ynms,zwsave)
          call sphtrans_fwd_cmpl_rsv(nterms,gx2mpoles(:,:,n,m),
     $    nphi,ntheta,gxs2d2grids(:,:,n,m),ctheta,whts,ynms,zwsave)
          call sphtrans_fwd_cmpl_rsv(nterms,gy2mpoles(:,:,n,m),
     $    nphi,ntheta,gys2d2grids(:,:,n,m),ctheta,whts,ynms,zwsave)
          call sphtrans_fwd_cmpl_rsv(nterms,gz2mpoles(:,:,n,m),
     $    nphi,ntheta,gzs2d2grids(:,:,n,m),ctheta,whts,ynms,zwsave)
        enddo
      enddo

cc    prepare upsampled x,y,zigrid, wsigrid for near singular integral
       iup = 4
       iterms = iup * nterms
       call grid_dim_init(iterms,itheta,iphi)

       allocate(x1igrid(iphi,itheta), y1igrid(iphi,itheta), 
     $ z1igrid(iphi,itheta), ws1igrid(iphi,itheta),h1igrid(iphi,itheta),
     $ n1igrid(iphi,itheta,3))
       allocate(x2igrid(iphi,itheta), y2igrid(iphi,itheta), 
     $ z2igrid(iphi,itheta), ws2igrid(iphi,itheta),h2igrid(iphi,itheta),
     $ n2igrid(iphi,itheta,3))
       allocate(sphfigrid(iphi,itheta,0:nterms,-nterms:nterms))
       allocate(whts_i(itheta))

      call nearSing_upsamp_prep(iup, iterms, itheta, iphi,
     $  x1mpole, y1mpole, z1mpole,
     $  x1igrid, y1igrid, z1igrid, ws1igrid, h1igrid, n1igrid)

      call nearSing_upsamp_prep(iup, iterms, itheta, iphi,
     $  x2mpole, y2mpole, z2mpole,
     $  x2igrid, y2igrid, z2igrid, ws2igrid, h2igrid, n2igrid)

      !!! load sphfigrid and whts_i from mod_upsamp
      do n = 0, nterms
        do m = -n, n
          do i = 1, itheta
            do j = 1, iphi
              sphfigrid(j,i,n,m) = up_sphfgrid(j,i,n,m,iup)
            enddo  
          enddo
        enddo
      enddo

      do i = 1, itheta
        whts_i(i) = up_whts(i,iup)
      enddo


      h_far = 5 * pi/(nterms + 1)
      allocate(sphf_s(0:nterms, -nterms:nterms))

      write(*,*) "start: regular quadr"
cccccc  regular quadrature over drop 2 for source on drop1 
      do k = 1,ntheta
        do kk = 1,nphi
          call close_dist(x1grid(kk,k),y1grid(kk,k),z1grid(kk,k),
     $          x2grid, y2grid, z2grid, x2mpole, y2mpole, z2mpole,
     $          dist_cl, sphf_s, xs, ys, zs) !! x,y,zs are the closest point on surface
          write(9410,*) "dist_cl, xs, ys, zs"
          write(9410,*) k, kk, dist_cl, xs, ys, zs

          do n = 0,nterms
            do m = -n,n
          !!write(*,*) "using rotating sph", n, m

              !!! integral for Y_nm 
              !! if (x1,y1,z1) is NOT in the close distance
              if (dist_cl .ge. h_far) then 

                call integral_eval_efield_cmpl(
     $          nterms,ntheta,nphi,
     $          x2grid,y2grid,z2grid,ws2grid,sphfgrid(:,:,n,m),whts,   !! regular quadr weight
     $          x1grid(kk,k),y1grid(kk,k),z1grid(kk,k), 
     $          n1grid(kk,k,:), cres)  

              else 
                !! call near_singlar
                write(*,*) "s1d2 k,kk,n,m",k,kk,n,m,dist_cl
                call nearSing_integral_efield_cmpl(
     $          iterms, itheta, iphi, nintpol, 
     $          x2igrid, y2igrid, z2igrid, ws2igrid, 
     $          sphfigrid(:,:,n,m), whts_i,
     $         gx2mpoles(:,:,n,m),gy2mpoles(:,:,n,m),gz2mpoles(:,:,n,m),
     $          dist_cl, sphf_s,  !! interpolate singluar integral
     $          xs, ys, zs, x1grid(kk,k),y1grid(kk,k),z1grid(kk,k),
     $          n1grid(kk,k,:), cres)

             endif
             cls1d2grids(kk,k,n,m) = cres
             write(9412,*) kk,k,n,m,cls1d2grids(kk,k,n,m)
            enddo
          enddo

        enddo
      enddo     
ccccccccccccccc


cccccc  regular quadrature over drop 1 for source on drop2 
      do k = 1,ntheta
        do kk = 1,nphi
          call close_dist(x2grid(kk,k),y2grid(kk,k),z2grid(kk,k),
     $          x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole,
     $          dist_cl, sphf_s, xs, ys, zs) !! x,y,zs are the closest point on surface
          write(9420,*) "dist_cl, xs, ys, zs"
          write(9420,*) k, kk, dist_cl, xs, ys, zs

          do n = 0,nterms
            do m = -n,n

              if (dist_cl .ge. h_far) then 

                call integral_eval_efield_cmpl(
     $          nterms,ntheta,nphi,
     $          x1grid,y1grid,z1grid,ws1grid,sphfgrid(:,:,n,m),whts,   !! regular quadr weight
     $          x2grid(kk,k),y2grid(kk,k),z2grid(kk,k), 
     $          n2grid(kk,k,:), cres)  

              else 
                !! call near_singlar
                write(*,*) "s2d1 k,kk,n,m",k,kk,n,m,dist_cl
                call nearSing_integral_efield_cmpl(
     $          iterms, itheta, iphi, nintpol, 
     $          x1igrid, y1igrid, z1igrid, ws1igrid, 
     $          sphfigrid(:,:,n,m), whts_i,
     $         gx1mpoles(:,:,n,m),gy1mpoles(:,:,n,m),gz1mpoles(:,:,n,m),
     $          dist_cl, sphf_s,  !! interpolate singluar integral
     $          xs, ys, zs, x2grid(kk,k),y2grid(kk,k),z2grid(kk,k),
     $          n2grid(kk,k,:), cres)

                write(9421,*) dist_cl, cres

             endif
             cls2d1grids(kk,k,n,m) = cres
             write(9421,*) kk,k,n,m,cls2d1grids(kk,k,n,m),
     $ cls2d1grids(kk,k,n,m) - cls1d2grids(kk,ntheta+1-k,n,m)
    
            enddo
          enddo
        enddo
      enddo     

ccccccccc Inner Product ccccccccccccc
      allocate(
     $ cls1d1mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms),
     $ cls1d2mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms),
     $ cls2d1mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms),
     $ cls2d2mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms))

      do n = 0,nterms
        do m = -n,n
           call sphtrans_fwd_cmpl_rsv(nterms,cls1d1mpoles(:,:,n,m),
     $    nphi,ntheta,cls1d1grids(:,:,n,m),ctheta,whts,ynms,zwsave)
           call sphtrans_fwd_cmpl_rsv(nterms,cls2d2mpoles(:,:,n,m),
     $    nphi,ntheta,cls2d2grids(:,:,n,m),ctheta,whts,ynms,zwsave)
           call sphtrans_fwd_cmpl_rsv(nterms,cls1d2mpoles(:,:,n,m),
     $    nphi,ntheta,cls1d2grids(:,:,n,m),ctheta,whts,ynms,zwsave)
           call sphtrans_fwd_cmpl_rsv(nterms,cls2d1mpoles(:,:,n,m),
     $    nphi,ntheta,cls2d1grids(:,:,n,m),ctheta,whts,ynms,zwsave)
        enddo
      enddo  

      allocate(Einfn1mpole(0:nterms,-nterms:nterms),
     $         Einfn2mpole(0:nterms,-nterms:nterms))
      call sphtrans_fwd_cmpl_rsv(nterms,Einfn1mpole(:,:),
     $    nphi,ntheta,cEinfn1grid(:,:),ctheta,whts,ynms,zwsave)
      call sphtrans_fwd_cmpl_rsv(nterms,Einfn2mpole(:,:),
     $    nphi,ntheta,cEinfn2grid(:,:),ctheta,whts,ynms,zwsave)

ccccccccccccccccc build linear system ccccccccccccccccccc

      allocate(cefield_mat(N2,N2))
      j = 1
      do n = 1,nterms !! start 0
        do m = -n, n
          i = 1
          do n1 = 1,nterms !! start 0
            do m1 = -n1, n1
              cefield_mat(i, j)=cls1d1mpoles(n1,m1,n,m)
              cefield_mat(i, j + N2/2)=cls1d2mpoles(n1,m1,n,m)
              cefield_mat(i + N2/2, j)=cls2d1mpoles(n1,m1,n,m)
              cefield_mat(i + N2/2, j + N2/2)=cls2d2mpoles(n1,m1,n,m) 
           ! if (i == 1) then
        write(9401,*) i, j, cls1d1mpoles(n1,m1,n,m)!, 
!     $ cls1d2mpoles(n1,m1,n,m),cls2d1mpoles(n1,m1,n,m),
!     $ cls2d2mpoles(n1,m1,n,m)
       ! endif
    
              i = i + 1    
            enddo
          enddo
          j = j + 1
        enddo
      enddo    

      do i = 1, N2
        !! note: factor 4*pi is because the spherical harmonic functions we used are not normalized
         cefield_mat(i,i) = cefield_mat(i,i)+dcmplx(0.5d0,0d0)  
!         write(15, *) 'i=',i,cap_lpl(i,i)
      enddo
      
      allocate(cEinf_col(N2))
      write(14,*) 'Y proj of E_inf*ngrid, ALL B terms in A * X = B'
      write(21,*) "only nonzero B terms, in A * X = B"
      i = 1
      do n = 1, nterms !! start 0
          do m = -n, n              
              cEinf_col(i) = Einfn1mpole(n, m)
              cEinf_col(i + N2/2) = Einfn2mpole(n, m)
             ! write(14,*) 'n',n,'m',m,'i=',i,cap_Einf(i)
              i = i + 1
          enddo
      enddo

      do i = 1, N2
        if (abs(cEinf_col(i)*dconjg(cEinf_col(i))).lt.1d-25) then
          cEinf_col(i) = 0d0
        else
          write(21,*)  'i=',i, cEinf_col(i)
        endif
      enddo
      
cccccccccccc
       write(19,*)  " all A terms, in A * X = B"   
       write(20,*)  " only nonzero A terms, in A * X = B"  
       do k = 1, N2
         do l = 1, N2
           !write(19,*) 'i=',k,'j=',l, cap_lpl(k,l)
           tmp = abs(cefield_mat(k,l)*dconjg(cefield_mat(k,l)))
           if (tmp.lt.1d-25) then 
             cefield_mat(k,l) = 0d0
           else
            write(20,*) 'i=',k,'j=',l, cefield_mat(k,l)
           endif
         enddo
       enddo

      !!!! special case for the first step sphere shape
      !! avoid all zero entry on the first line of linear system (A and B) 
       if (istep .eq. 1) then 
        !cefield_mat(1,1) = (1d0, 0d0)
        !cefield_mat(1+N2/2,1+N2/2) = (1d0, 0d0)
      endif                 
       

cccccccccccc
      allocate(ipiv(N2))         !! n start with 0            
      CALL zgesv(N2,NRHS,cefield_mat,N2,IPIV,cEinf_col,N2,INFO)      
!      CALL zgesv(N2-1,NRHS,cap_lpl,N2-1,
!     $ IPIV,cap_Einf,N2-1,INFO)      
      write(16,*) 'linear solver zgesv solution:','INFO=',INFO

      i = 1
      do n = 1, nterms !! start 0
          do m = -n, n              
              en1mpole(n, m) = cEinf_col(i)   
              en2mpole(n, m) = cEinf_col(i + N2/2)
              write(16,*) 'i=',i,'n',n,'m',m,cEinf_col(i)
              tmp = abs(en1mpole(n, m)*dconjg(en1mpole(n, m)))
              if (tmp.lt.1d-25) then
                en1mpole(n, m) = 0d0
              endif
              tmp = abs(en2mpole(n, m)*dconjg(en2mpole(n, m)))
              if (tmp.lt.1d-25) then
                en2mpole(n, m) = 0d0
              endif
              i = i + 1
          enddo
      enddo

      en1mpole(0,0) = (0d0, 0d0)
      en2mpole(0,0) = (0d0, 0d0)

!      !!!! special case for the first step sphere shape
!      if (istep .eq. 1) then 
!        enmpole(0,0) = (0d0, 0d0)
!      endif
      write(111,*) "n, m, en1mpole(n,m), en2mpole(n,m)"
      do n = 1, nterms
        do m = -n, n
          write(111,*) n, m, en1mpole(n,m), en2mpole(n,m)
          do i = 1, ntheta
            do j = 1, nphi
          write(112,*) i,j,n,m,sphfgrid(j,i,n,m)
            enddo
          enddo
        enddo
      enddo      
ccccc spherical harmonic transform ccccccc
      allocate(En1(nphi, ntheta),En2(nphi, ntheta))

      do i = 1,ntheta
          do j = 1, nphi
            ctmp1 = (0d0,0d0)
            ctmp2 = (0d0,0d0)
            do n = 0, nterms   !! start 0, only take odd terms
              do m = -n, n              
                  ctmp1 = ctmp1 + en1mpole(n, m)*sphfgrid(j,i,n,m) 
                  ctmp2 = ctmp2 + en2mpole(n, m)*sphfgrid(j,i,n,m)                
              enddo
            enddo 
            En1(j,i) = dble(ctmp1)
            En2(j,i) = dble(ctmp2)
            write(17,*) 'j=',j,'i=',i, En1(j,i),En2(j,i)
          enddo
      enddo
cccccccccc total force mpole cccccccccc
c          ( En^2/2*Ca + 2*H ) * ngrid
ccccccccccccccccccccccccccccccccccccccc    
       write(*,*) "start: evaluate velocity"

      allocate(fx1(nphi,ntheta), fy1(nphi,ntheta), fz1(nphi,ntheta),
     $ fx2(nphi,ntheta), fy2(nphi,ntheta), fz2(nphi,ntheta))

      Caflow = 0.0d0
      do i = 1, ntheta
        do j = 1, nphi            
          tmp1 = (En1(j,i)**2/2d0*Ca + 2d0*h1grid(j,i))
          fx1(j,i) = tmp1*n1grid(j,i,1)
          fy1(j,i) = tmp1*n1grid(j,i,2)
          fz1(j,i) = tmp1*n1grid(j,i,3)
          tmp2 = (En2(j,i)**2/2d0*Ca + 2d0*h2grid(j,i))
          fx2(j,i) = tmp2*n2grid(j,i,1)
          fy2(j,i) = tmp2*n2grid(j,i,2)
          fz2(j,i) = tmp2*n2grid(j,i,3)
          !write(32,*) j,i,En(j,i)/dcos(theta(i)),dcos(theta(i)),ws(j,i)
        enddo
      enddo

      allocate(fx1mpole(0:nterms,0:nterms), fy1mpole(0:nterms,0:nterms),
     $ fz1mpole(0:nterms,0:nterms),fx2mpole(0:nterms,0:nterms),
     $ fy2mpole(0:nterms,0:nterms),fz2mpole(0:nterms,0:nterms))
      
      call sphtrans_fwd_real_rsv(nterms,fx1mpole,nphi,ntheta,fx1,
     $     ctheta,whts,ynms,dwsave)
      call sphtrans_fwd_real_rsv(nterms,fy1mpole,nphi,ntheta,fy1,
     $     ctheta,whts,ynms,dwsave)    
      call sphtrans_fwd_real_rsv(nterms,fz1mpole,nphi,ntheta,fz1,
     $     ctheta,whts,ynms,dwsave)           
      call sphtrans_fwd_real_rsv(nterms,fx2mpole,nphi,ntheta,fx2,
     $     ctheta,whts,ynms,dwsave)
      call sphtrans_fwd_real_rsv(nterms,fy2mpole,nphi,ntheta,fy2,
     $     ctheta,whts,ynms,dwsave)    
      call sphtrans_fwd_real_rsv(nterms,fz2mpole,nphi,ntheta,fz2,
     $     ctheta,whts,ynms,dwsave)           
cccccccccccccccccccccccccccccccccccccccc

cccccc rotate force components mpole ccccccccccc

      allocate(rfx1grids(nphi,ntheta,nrot,nbeta),
     $ rfy1grids(nphi,ntheta,nrot,nbeta),
     $ rfz1grids(nphi,ntheta,nrot,nbeta),
     $ rfx2grids(nphi,ntheta,nrot,nbeta),
     $ rfy2grids(nphi,ntheta,nrot,nbeta),
     $ rfz2grids(nphi,ntheta,nrot,nbeta))

      call rotgrid_dsr_real(nterms,fx1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rfx1grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,fy1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rfy1grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,fz1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rfz1grids,
     $ rotmat,ctheta,ynms,dwsave)            
      call rotgrid_dsr_real(nterms,fx2mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rfx2grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,fy2mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rfy2grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,fz2mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rfz2grids,
     $ rotmat,ctheta,ynms,dwsave)    
cccccccccccccccccccccccccccccccccccccccccccccc
      
      allocate(us1d1grid(nphi,ntheta,3),us2d2grid(nphi,ntheta,3),
     $ us1d2grid(nphi,ntheta,3),us2d1grid(nphi,ntheta,3))
cc    special gauss-legendre quadrature(lapalce) for all rotated grids
          write(34,*) 'irot ','jbeta ', 'special g-l quadr stk'
          do kk = 1,nrot
            do k = 1,nbeta
              call integral_eval_ufield(ntheta,nphi,
     $ rx1grids(:,:,kk,k), ry1grids(:,:,kk,k), rz1grids(:,:,kk,k), 
     $ rws1grids(:,:,kk,k), rfx1grids(:,:,kk,k), rfy1grids(:,:,kk,k), 
     $ rfz1grids(:,:,kk,k), swhts, 
     $ x1grid(kk,k), y1grid(kk,k), z1grid(kk,k), 
     $ quadr1, quadr2, quadr3)
              us1d1grid(kk,k,1) = quadr1 
              us1d1grid(kk,k,2) = quadr2
              us1d1grid(kk,k,3) = quadr3 
              !write(34,*) kk, k, us1d1grid(kk,k,1), us1d1grid(kk,k,2),
   !  $  us1d1grid(kk,k,3)
            enddo
          enddo 

          do kk = 1,nrot
            do k = 1,nbeta
              call integral_eval_ufield(ntheta,nphi,
     $ rx2grids(:,:,kk,k), ry2grids(:,:,kk,k), rz2grids(:,:,kk,k), 
     $ rws2grids(:,:,kk,k), rfx2grids(:,:,kk,k), rfy2grids(:,:,kk,k), 
     $ rfz2grids(:,:,kk,k), swhts, 
     $ x2grid(kk,k), y2grid(kk,k), z2grid(kk,k), 
     $ quadr1, quadr2, quadr3)
              us2d2grid(kk,k,1) = quadr1
              us2d2grid(kk,k,2) = quadr2 
              us2d2grid(kk,k,3) = quadr3
             ! write(34,*) kk, k, us2d2grid(kk,k,1), us2d2grid(kk,k,2),
    ! $  us2d2grid(kk,k,3)
            enddo
          enddo 

      allocate(u1mpoles(0:nterms,0:nterms,3),
     $ u2mpoles(0:nterms,0:nterms,3))

      do i = 1, 3
             call sphtrans_fwd_real_rsv(nterms,u1mpoles(:,:,i),
     $    nphi,ntheta,us1d1grid(:,:,i),ctheta,whts,ynms,dwsave)
             call sphtrans_fwd_real_rsv(nterms,u2mpoles(:,:,i),
     $    nphi,ntheta,us2d2grid(:,:,i),ctheta,whts,ynms,dwsave)      
      enddo


      allocate(fx1igrid(iphi,itheta), fy1igrid(iphi,itheta),
     $ fz1igrid(iphi,itheta), fx2igrid(iphi,itheta), 
     $ fy2igrid(iphi,itheta), fz2igrid(iphi,itheta))

      allocate(En1igrid(iphi,itheta), En2igrid(iphi,itheta))

      call nearSing_upsamp_En(iup, iterms, itheta, iphi,
     $  en1mpole, en2mpole, En1igrid, En2igrid)
      
      do i = 1, itheta
        do j = 1, iphi            
          tmp1 = (En1igrid(j,i)**2/2d0*Ca + 2d0*h1igrid(j,i))
          fx1igrid(j,i) = tmp1*n1igrid(j,i,1)
          fy1igrid(j,i) = tmp1*n1igrid(j,i,2)
          fz1igrid(j,i) = tmp1*n1igrid(j,i,3)
          tmp2 = (En2igrid(j,i)**2/2d0*Ca + 2d0*h2igrid(j,i))
          fx2igrid(j,i) = tmp2*n2igrid(j,i,1)
          fy2igrid(j,i) = tmp2*n2igrid(j,i,2)
          fz2igrid(j,i) = tmp2*n2igrid(j,i,3)
          !write(32,*) j,i,En(j,i)/dcos(theta(i)),dcos(theta(i)),ws(j,i)
        enddo
      enddo


          do kk = 1,nphi
            do k = 1,ntheta
              !!! integral for Y_nm 
              !! if (x1,y1,z1) is NOT in the close distance
                 
                call close_dist(x1grid(kk,k),y1grid(kk,k),z1grid(kk,k),
     $          x2grid, y2grid, z2grid, x2mpole, y2mpole, z2mpole,
     $          dist_cl, sphf_s, xs, ys, zs) !! x,y,zs are the closest point on surface
  
              if (dist_cl .ge. h_far) then 

                call integral_eval_ufield(ntheta, nphi,
     $          x2grid, y2grid, z2grid, 
     $          ws2grid, fx2, fy2, fz2, whts,     !! regular whts
     $          x1grid(kk,k), y1grid(kk,k), z1grid(kk,k), 
     $          quadr1, quadr2, quadr3)

              else 
                !! call near_singlar
                call nearSing_integral_ufield(
     $          iterms, itheta, iphi, nintpol, 
     $          x2igrid, y2igrid, z2igrid, ws2igrid, 
     $          fx2igrid, fy2igrid, fz2igrid, whts_i,
     $          u2mpoles, dist_cl, sphf_s,  !! interpolate singluar integral
     $          xs, ys, zs, x1grid(kk,k),y1grid(kk,k),z1grid(kk,k),
     $          quadr1, quadr2, quadr3)

             endif
             us1d2grid(kk,k,1) = quadr1
             us1d2grid(kk,k,2) = quadr2 
             us1d2grid(kk,k,3) = quadr3
            enddo
          enddo


          do kk = 1,nphi
            do k = 1,ntheta
              !!! integral for Y_nm 
              !! if (x1,y1,z1) is NOT in the close distance
                 
                call close_dist(x2grid(kk,k),y2grid(kk,k),z2grid(kk,k),
     $          x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole,
     $          dist_cl, sphf_s, xs, ys, zs) !! x,y,zs are the closest point on surface
  
              if (dist_cl .ge. h_far) then 

                call integral_eval_ufield(ntheta, nphi,
     $          x1grid, y1grid, z1grid, 
     $          ws1grid, fx1, fy1, fz1, whts,     !! regular whts
     $          x2grid(kk,k), y2grid(kk,k), z2grid(kk,k), 
     $          quadr1, quadr2, quadr3)

              else 
                !! call near_singlar
                call nearSing_integral_ufield(
     $          iterms, itheta, iphi, nintpol, 
     $          x1igrid, y1igrid, z1igrid, ws1igrid, 
     $          fx1igrid, fy1igrid, fz1igrid, whts_i,
     $          u1mpoles, dist_cl, sphf_s,  !! interpolate singluar integral
     $          xs, ys, zs, x2grid(kk,k), y2grid(kk,k), z2grid(kk,k),
     $          quadr1, quadr2, quadr3)

             endif
             us2d1grid(kk,k,1) = quadr1
             us2d1grid(kk,k,2) = quadr2 
             us2d1grid(kk,k,3) = quadr3
            enddo
          enddo    
          
          do i = 1, ntheta
            do j = 1, nphi
              do k = 1,3
              u1grid(j,i,k) = us1d1grid(j,i,k) + us1d2grid(j,i,k)
              u2grid(j,i,k) = us2d1grid(j,i,k) + us2d2grid(j,i,k)
              enddo
            enddo
          enddo

!!!<       electrical total force on drop 1
          ele_tot_fx = 0d0
          ele_tot_fy = 0d0
          ele_tot_fz = 0d0
          do i = 1, ntheta
            tmpx = 0d0
            tmpy = 0d0
            tmpz = 0d0
            do j = 1, nphi            
              tmp = En1(j,i)**2/2d0!*Ca
              tmpx = tmpx + tmp*n1grid(j,i,1)
              tmpy = tmpy + tmp*n1grid(j,i,2)
              tmpz = tmpz + tmp*n1grid(j,i,3)
            enddo
            ele_tot_fx = ele_tot_fx + tmpx * whts(i)
            ele_tot_fy = ele_tot_fy + tmpy * whts(i)
            ele_tot_fz = ele_tot_fz + tmpz * whts(i)
          enddo
          ele_tot_fx = ele_tot_fx/nphi*(2*pi)
          ele_tot_fy = ele_tot_fy/nphi*(2*pi)
          ele_tot_fz = ele_tot_fz/nphi*(2*pi)
          ele_tot_fmag = ele_tot_fx**2 + ele_tot_fy**2 + ele_tot_fz**2
          ele_tot_fmag = dsqrt(ele_tot_fmag)
          write(109,*) "electrical total force on drop 1"
          write(109,*) "x,y,z component and magnitude"
          write(109,*) ele_tot_fx, ele_tot_fy, ele_tot_fz, ele_tot_fmag
!!!>



      deallocate(x1igrid, y1igrid, z1igrid, ws1igrid, h1igrid,
     $ x2igrid, y2igrid, z2igrid, ws2igrid, h2igrid, whts_i)
      deallocate(n1igrid, n2igrid, sphfigrid, cEinfn1grid,cEinfn2grid)
      deallocate(ws1grid, ws2grid, ws1mpole,ws2mpole,rx1grids, ry1grids,
     $ rz1grids, rws1grids, rx2grids, ry2grids, rz2grids, rws2grids)
     !!deallocate(f1mpoles,f2mpoles)
      deallocate(cls1d1grids, cls1d2grids, cls2d1grids, cls2d2grids)
      deallocate(cls1d1mpoles,cls1d2mpoles, cls2d1mpoles, cls2d2mpoles)
      deallocate(Einfn1mpole, Einfn2mpole, cefield_mat, cEinf_col, ipiv)
      deallocate(En1, En2, fx1, fy1, fz1, fx2, fy2, fz2)
      deallocate(fx1mpole, fy1mpole, fz1mpole, fx2mpole, fy2mpole, 
     $ fz2mpole)
      deallocate(rfx1grids, rfy1grids, rfz1grids, rfx2grids, rfy2grids, 
     $ rfz2grids)
      deallocate(us1d1grid, us2d1grid, us1d2grid, us2d2grid)
      deallocate(u1mpoles, u2mpoles)
      deallocate(fx1igrid, fy1igrid, fz1igrid, fx2igrid, fy2igrid, 
     $ fz2igrid, En1igrid, En2igrid)

      end SUBROUTINE

c
c  Our definition of complex spherical harmonics is
c
c  Ynm(theta,phi)= sqrt(2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(im phi), 
c  Yn,-m(theta,phi) = sqrt(2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(-im phi),   for m >= 0.
c       
c  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.




 
