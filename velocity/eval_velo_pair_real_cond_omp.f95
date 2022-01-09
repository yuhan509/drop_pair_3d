SUBROUTINE eval_velo_pair_cond_wmod(istep,  &
    x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole, h1mpole, ws1mpole, &
    n1mpole, h1grid, ws1grid, n1grid, u1grid, en1mpole,   &
    x2grid, y2grid, z2grid, x2mpole, y2mpole, z2mpole, h2mpole, ws2mpole, &
    n2mpole, h2grid, ws2grid, n2grid, u2grid, en2mpole)
  
  use OMP_LIB
  use mod_orgHelper
  use mod_rotsphf
  use mod_upsamp
       
  implicit real *8 (a-h,o-z)
       
  integer :: istat = 2, istep, iup
  
  real *8 :: x1grid(nphi,ntheta),y1grid(nphi,ntheta),z1grid(nphi,ntheta), &
  h1grid(nphi,ntheta),ws1grid(nphi,ntheta), n1grid(nphi,ntheta,3)
  complex *16 :: x1mpole(0:nterms,0:nterms),y1mpole(0:nterms,0:nterms),z1mpole(0:nterms,0:nterms)
  complex *16 :: h1mpole(0:nterms,0:nterms),ws1mpole(0:nterms,0:nterms),n1mpole(0:nterms,0:nterms,3)
  
  real *8 :: x2grid(nphi,ntheta),y2grid(nphi,ntheta),z2grid(nphi,ntheta), &
  h2grid(nphi,ntheta),ws2grid(nphi,ntheta), n2grid(nphi,ntheta,3)
  complex *16 :: x2mpole(0:nterms,0:nterms),y2mpole(0:nterms,0:nterms),z2mpole(0:nterms,0:nterms)
  complex *16 :: h2mpole(0:nterms,0:nterms),ws2mpole(0:nterms,0:nterms),n2mpole(0:nterms,0:nterms,3)
  
       !! output 
  real *8 :: u1grid(nphi,ntheta,3),u2grid(nphi,ntheta,3)
  complex *16 :: en1mpole(0:nterms,0:nterms),en2mpole(0:nterms,0:nterms)
       !! end of output
       
       complex *16, ALLOCATABLE :: sphf_s(:,:)
       integer :: iterms, itheta, iphi
       real *8, ALLOCATABLE :: whts_i(:)
       real *8, dimension(:,:),ALLOCATABLE ::  &
     x1igrid, y1igrid, z1igrid, ws1igrid, h1igrid,    &
     x2igrid, y2igrid, z2igrid, ws2igrid, h2igrid
       real *8, dimension(:,:,:),ALLOCATABLE :: n1igrid, n2igrid
  
       complex *16,dimension(:,:),ALLOCATABLE:: cpinf1grid, cpinf2grid
     !  real *8, dimension(:,:), ALLOCATABLE:: ws1grid, ws2grid
     !  complex *16, dimension(:,:), ALLOCATABLE:: ws1mpole, ws2mpole
       real *8, dimension(:,:,:,:), ALLOCATABLE:: rx1grids, ry1grids, &
      rz1grids, rws1grids, rx2grids, ry2grids, rz2grids, rws2grids
  
       complex *16,dimension(:,:,:,:),ALLOCATABLE :: g1mpoles,g2mpoles
  
       complex *16,dimension(:,:,:,:),ALLOCATABLE :: f1mpoles,f2mpoles
  
       complex *16, dimension(:,:,:,:), ALLOCATABLE :: cls1d1grids, &
  cls1d2grids, cls2d1grids, cls2d2grids
       complex *16, dimension(:,:,:,:), ALLOCATABLE :: cls1d1mpoles, &
  cls1d2mpoles, cls2d1mpoles, cls2d2mpoles
  
       complex *16, dimension(:,:), ALLOCATABLE:: Pinf1mpole,Pinf2mpole
       complex *16, dimension(:,:), ALLOCATABLE:: cdrop1_chrg_LHS,cdrop2_chrg_LHS

       complex *16, ALLOCATABLE:: cLHS_mat(:,:), cRHS_col(:) 
  
       INTEGER :: NRHS=1,N2,N21
       INTEGER :: INFO
       INTEGER, allocatable :: ipiv(:)
  
       real *8, ALLOCATABLE, dimension(:,:):: En1, En2
       real *8, dimension(:,:), allocatable:: fx1,fy1,fz1,fx2,fy2,fz2
       complex *16, dimension(:,:), allocatable :: fx1mpole, fy1mpole, &
   fz1mpole, fx2mpole, fy2mpole, fz2mpole
       real *8, dimension(:,:,:,:), allocatable :: rfx1grids, rfy1grids, &
   rfz1grids, rfx2grids, rfy2grids, rfz2grids
  
       real *8, dimension(:,:,:), allocatable :: us1d1grid, us2d1grid,  & 
   us1d2grid, us2d2grid
       complex *16, dimension(:,:,:), allocatable :: u1mpoles, u2mpoles
       real *8, dimension(:,:), allocatable :: fx1igrid, fy1igrid, & 
   fz1igrid, fx2igrid, fy2igrid, fz2igrid, En1igrid, En2igrid
       
       complex* 16::ctmp,cquadr,ctmp1,ctmp2,ctmp3,cquadr1,cquadr2,cquadr3,cres
       DOUBLE PRECISION START, END  !! time stamp
       !! dimension of matrix block (exclude n=0,m=0 in expansions for En1, En2)
       N2 = ((nterms+1)*(nterms+1))*2 + 2
       N21 = (nterms+1)*(nterms+1)
  
       PI = 4.D0*DATAN(1.D0)
       iup = maxup   !! maxup is assigned in mod_assign
  
  !cccccccccccccccccccccccccccccccccccccccccccccc
       !allocate(ws1grid(nphi,ntheta), ws2grid(nphi,ntheta))
       allocate(cpinf1grid(nphi,ntheta),cpinf2grid(nphi,ntheta))
       
       do i=1,ntheta
         do j=1,nphi
          ! ws1grid(j,i) = w1grid(j,i)/dsin(theta(i))
          ! ws2grid(j,i) = w2grid(j,i)/dsin(theta(i))
           !! E_inf * ngrid
           cpinf1grid(j,i) = dcmplx(z1grid(j,i)*dcos(alpha) + &
               x1grid(j,i)*dsin(alpha), 0d0)
           cpinf2grid(j,i) = dcmplx(z2grid(j,i)*dcos(alpha) + &
               x2grid(j,i)*dsin(alpha), 0d0)
           !write(88,*) j, i, ngrid(j,i,3), ctheta(i)
         enddo
       enddo
       write(*,*) "Pinfn assign done"
  
   !    allocate(ws1mpole(0:nterms,0:nterms),ws2mpole(0:nterms,0:nterms))
  
   !    call sphtrans_fwd_real_rsv(nterms,ws1mpole,nphi,ntheta,ws1grid,
   !   $  ctheta,whts,ynms,dwsave)
   !    call sphtrans_fwd_real_rsv(nterms,ws2mpole,nphi,ntheta,ws2grid,
   !   $  ctheta,whts,ynms,dwsave)
  !cccccccccccccccccccccccccccccccccccc
  !cccccc rotate x,y,z mpole and local area element ccccccccccc
       nbeta = ntheta
       nrot = nphi
  
  allocate(rx1grids(nphi,ntheta,nrot,nbeta), ry1grids(nphi,ntheta,nrot,nbeta), &
   rz1grids(nphi,ntheta,nrot,nbeta), rws1grids(nphi,ntheta,nrot,nbeta))
  
  allocate(rx2grids(nphi,ntheta,nrot,nbeta), ry2grids(nphi,ntheta,nrot,nbeta), &
   rz2grids(nphi,ntheta,nrot,nbeta), rws2grids(nphi,ntheta,nrot,nbeta))
  
  call rotgrid_fsr_real(nterms,x1mpole,nphi,ntheta,nbeta,theta,nrot,rx1grids,rotmat,ctheta,ynms,dwsave)
  call rotgrid_fsr_real(nterms,y1mpole,nphi,ntheta,nbeta,theta,nrot,ry1grids,rotmat,ctheta,ynms,dwsave)
  call rotgrid_fsr_real(nterms,z1mpole,nphi,ntheta,nbeta,theta,nrot,rz1grids,rotmat,ctheta,ynms,dwsave)     
  call rotgrid_fsr_real(nterms,ws1mpole,nphi,ntheta,nbeta,theta,nrot,rws1grids,rotmat,ctheta,ynms,dwsave)     
  
  call rotgrid_fsr_real(nterms,x2mpole,nphi,ntheta,nbeta,theta,nrot,rx2grids,rotmat,ctheta,ynms,dwsave)
  call rotgrid_fsr_real(nterms,y2mpole,nphi,ntheta,nbeta,theta,nrot,ry2grids,rotmat,ctheta,ynms,dwsave)
  call rotgrid_fsr_real(nterms,z2mpole,nphi,ntheta,nbeta,theta,nrot,rz2grids,rotmat,ctheta,ynms,dwsave)     
  call rotgrid_fsr_real(nterms,ws2mpole,nphi,ntheta,nbeta,theta,nrot,rws2grids,rotmat,ctheta,ynms,dwsave)     
  
  
       allocate(cls1d1grids(nrot,nbeta,0:nterms,-nterms:nterms),&
         cls2d2grids(nrot,nbeta,0:nterms,-nterms:nterms), &
         cls1d2grids(nrot,nbeta,0:nterms,-nterms:nterms), &
         cls2d1grids(nrot,nbeta,0:nterms,-nterms:nterms))
  
       write(*,*) "start: special quadr efield"
  
  !cc    special gauss-legendre quadrature(lapalce) for all rotated grids
       !$OMP PARALLEL 
       !$OMP do schedule(dynamic) private(kk,n,m,cquadr1, cquadr2, cquadr3)  
       do k = 1,nbeta
         do kk = 1,nrot
           !!write(*,*) "using rotating sph", n, m
           do n = 0,nterms
             do m = 0,n  !!!real, change to zero
              ! write(*,*) "s1d1 n,m,kk,k", n,m,kk,k
               call integral_g_efield_cmpl(nterms,ntheta,nphi, &
                       rx1grids(:,:,kk,k),ry1grids(:,:,kk,k),rz1grids(:,:,kk,k), &
                       rws1grids(:,:,kk,k), rsphfgrids(:,:,kk,k,n,m), swhts, &
                       x1grid(kk,k), y1grid(kk,k), z1grid(kk,k),  &
                       cquadr1)  
               !! naming : c- is for complex; l- is for laplace kernel;
               !! s1- means source(kernel) point is on drop 1;
               !! d1- means the contribution from integrating over drop 1.
               !! (kk,k) represents kernel point, will serve as 
               cls1d1grids(kk,k,n,m) = cquadr1
              ! write(9411,*) kk,k,n,m,cls1d1grids(kk,k,n,m)
  
                  ! write(*,*) "s2d2 n,m,kk,k", n,m,kk,k
              call integral_g_efield_cmpl(nterms,ntheta,nphi, &
                      rx2grids(:,:,kk,k),ry2grids(:,:,kk,k),rz2grids(:,:,kk,k), &
                      rws2grids(:,:,kk,k), rsphfgrids(:,:,kk,k,n,m), swhts, &
                      x2grid(kk,k), y2grid(kk,k), z2grid(kk,k),  &
                      cquadr2)  
     
               cls2d2grids(kk,k,n,m) = cquadr2 
                  ! write(9422,*) kk,k,n,m,cls2d2grids(kk,k,n,m)
             enddo
           enddo
         enddo
       enddo     
       !$OMP end do
       !$OMP END PARALLEL 
  
  
       allocate(g1mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms), &
       g2mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms))
      
       do n = 0, nterms
         do m = 0, n !!! real, change to 0
           call sphtrans_fwd_cmpl_rsv(nterms,g1mpoles(:,:,n,m), &
              nphi,ntheta,cls1d1grids(:,:,n,m),ctheta,whts,ynms,zwsave)
           call sphtrans_fwd_cmpl_rsv(nterms,g2mpoles(:,:,n,m), &
              nphi,ntheta,cls2d2grids(:,:,n,m),ctheta,whts,ynms,zwsave)
         enddo
       enddo
  
  !cc    prepare upsampled x,y,zigrid, wsigrid for near singular integral
  
        iterms = iup * nterms
        call grid_dim_init(iterms,itheta,iphi)
  
        allocate(x1igrid(iphi,itheta), y1igrid(iphi,itheta), &
         z1igrid(iphi,itheta), ws1igrid(iphi,itheta),h1igrid(iphi,itheta),&
         n1igrid(iphi,itheta,3))
        allocate(x2igrid(iphi,itheta), y2igrid(iphi,itheta), &
         z2igrid(iphi,itheta), ws2igrid(iphi,itheta),h2igrid(iphi,itheta),&
         n2igrid(iphi,itheta,3))
        allocate(whts_i(itheta))
  
       call nearSing_upsamp_prep_mpole(iup, iterms, itheta, iphi,&
         x1mpole, y1mpole, z1mpole, h1mpole, ws1mpole, n1mpole, &
         x1igrid, y1igrid, z1igrid, ws1igrid, h1igrid, n1igrid)
  
       call nearSing_upsamp_prep_mpole(iup, iterms, itheta, iphi,&
         x2mpole, y2mpole, z2mpole, h2mpole, ws2mpole, n2mpole, &
         x2igrid, y2igrid, z2igrid, ws2igrid, h2igrid, n2igrid)
  
       do i = 1, itheta
         whts_i(i) = up_whts(i,iup)
        ! write(9405,*) up_whts(i,iup)
       enddo
  
       allocate(sphf_s(0:nterms, -nterms:nterms))
  
       write(*,*) "start: regular quadr"

       iflag_cd = 2
  !cccccc  regular quadrature over drop 2 for source on drop1 
       START = omp_get_wtime() 
       !$OMP PARALLEL 
       !$OMP do schedule(static) private(k,n,m,dist_cl, sphf_s, xs, ys, zs, cres)  
       do kk = 1,nphi
         do k = 1,ntheta
  
           call close_dist_flag(iflag_cd, x1grid(kk,k),y1grid(kk,k),z1grid(kk,k), &
                     x2grid, y2grid, z2grid, x2mpole, y2mpole, z2mpole, &
                     dist_cl, sphf_s, xs, ys, zs) !! x,y,zs are the closest point on surface
           !write(9410,*) "dist_cl, xs, ys, zs"
           !write(9410,*) k, kk, dist_cl, xs, ys, zs
  
           do n = 0,nterms
             do m = 0,n   !!! real, change to 0
           !!write(*,*) "using rotating sph", n, m
  
               !!! integral for Y_nm 
               !! if (x1,y1,z1) is NOT in the close distance
               if (dist_cl .ge. h_far) then 
  
                 call integral_g_efield_cmpl( &
                           nterms,ntheta,nphi, &
                           x2grid,y2grid,z2grid,ws2grid,sphfgrid(:,:,n,m),whts, &  !! regular quadr weight 
                           x1grid(kk,k),y1grid(kk,k),z1grid(kk,k), cres)  
  
               else 
                 !! call near_singlar
                 !write(*,*) "s1d2 k,kk,n,m",k,kk,n,m,dist_cl
                 call nearSing_integral_g_efield_cmpl( iterms, itheta, iphi,  &
                           x2igrid, y2igrid, z2igrid, ws2igrid,  &
                           upx_sphfgrid(:,:,n,m), whts_i, g2mpoles(:,:,n,m), &
                           dist_cl, sphf_s, & !! interpolate singluar integral 
                           xs, ys, zs, x1grid(kk,k),y1grid(kk,k),z1grid(kk,k), cres)
  
              endif
              cls1d2grids(kk,k,n,m) = cres

!              if (n .eq. 0) then
!              write(9910,*) "source: k,kk", k, kk, cres
!              endif
             enddo
           enddo
  
       ! source point on dorp2 integral over drop1
          call close_dist_flag(iflag_cd, x2grid(kk,k),y2grid(kk,k),z2grid(kk,k), &
                     x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole, &
                     dist_cl, sphf_s, xs, ys, zs) !! x,y,zs are the closest point on surface
           !write(9420,*) "dist_cl, xs, ys, zs"
           !write(9420,*) k, kk, dist_cl, xs, ys, zs
  
           do n = 0,nterms
             do m = 0,n !!! real, change to 0
  
               if (dist_cl .ge. h_far) then 
  
                 call integral_g_efield_cmpl( &
                           nterms,ntheta,nphi, &
                           x1grid,y1grid,z1grid,ws1grid,sphfgrid(:,:,n,m),whts, &   !! regular quadr weight
                           x2grid(kk,k),y2grid(kk,k),z2grid(kk,k), cres)  
  
               else 
                 !! call near_singlar
                ! write(*,*) "s2d1 k,kk,n,m",k,kk,n,m,dist_cl
                 call nearSing_integral_g_efield_cmpl(iterms, itheta, iphi, &
                           x1igrid, y1igrid, z1igrid, ws1igrid,  &
                           upx_sphfgrid(:,:,n,m), whts_i, g1mpoles(:,:,n,m), &
                           dist_cl, sphf_s, & !! interpolate singluar integral 
                           xs, ys, zs, x2grid(kk,k),y2grid(kk,k),z2grid(kk,k), cres)
  
                ! write(9421,*) dist_cl, cres
  
              endif
              !if (m .eq. 0) then
               cls2d1grids(kk,k,n,m) = cres
              ! else 
              ! cls2d1grids(kk,k,n,m) = cres * 2
              !endif
            !  write(9421,*) kk,k,n,m,cls2d1grids(kk,k,n,m),
     ! $ cls2d1grids(kk,k,n,m) - cls1d2grids(kk,ntheta+1-k,n,m)
     
             enddo
           enddo 
  
         enddo
       enddo     
       !$OMP end do
       !$OMP END PARALLEL 
  !ccccccccccccccc
       END = omp_get_wtime() 
       PRINT *, "Work took", END - START, "seconds"       
  
  !ccccccccc Inner Product ccccccccccccc
       allocate(cls1d1mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms), &
       cls1d2mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms), &
       cls2d1mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms), &
       cls2d2mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms))
  
       do n = 0,nterms
         do m = 0,n  !!! real, change to 0
            call sphtrans_fwd_cmpl_rsv(nterms,cls1d1mpoles(:,:,n,m), &
                nphi,ntheta,cls1d1grids(:,:,n,m),ctheta,whts,ynms,zwsave)
            call sphtrans_fwd_cmpl_rsv(nterms,cls2d2mpoles(:,:,n,m), &
                nphi,ntheta,cls2d2grids(:,:,n,m),ctheta,whts,ynms,zwsave)
            call sphtrans_fwd_cmpl_rsv(nterms,cls1d2mpoles(:,:,n,m), &
                nphi,ntheta,cls1d2grids(:,:,n,m),ctheta,whts,ynms,zwsave)
            call sphtrans_fwd_cmpl_rsv(nterms,cls2d1mpoles(:,:,n,m), &
                nphi,ntheta,cls2d1grids(:,:,n,m),ctheta,whts,ynms,zwsave)
         enddo
       enddo
       
       do n = 1, nterms
         do m = 1, n
           do n1 = 0, nterms
             do m1 = -n1, n1
              cls1d1mpoles(n1,m1,n,-m) = dconjg(cls1d1mpoles(n1,-m1,n,m))
              cls1d2mpoles(n1,m1,n,-m) = dconjg(cls1d2mpoles(n1,-m1,n,m))
              cls2d1mpoles(n1,m1,n,-m) = dconjg(cls2d1mpoles(n1,-m1,n,m))
              cls2d2mpoles(n1,m1,n,-m) = dconjg(cls2d2mpoles(n1,-m1,n,m))
             enddo
           enddo
         enddo
       enddo
  
       allocate(Pinf1mpole(0:nterms,-nterms:nterms), &
               Pinf2mpole(0:nterms,-nterms:nterms))

       call sphtrans_fwd_cmpl_rsv(nterms,Pinf1mpole(:,:), &
          nphi,ntheta,cpinf1grid(:,:),ctheta,whts,ynms,zwsave)
       call sphtrans_fwd_cmpl_rsv(nterms,Pinf2mpole(:,:), &
          nphi,ntheta,cpinf2grid(:,:),ctheta,whts,ynms,zwsave)
  
       allocate(cdrop1_chrg_LHS(0:nterms,-nterms:nterms), &
       cdrop2_chrg_LHS(0:nterms,-nterms:nterms))

          do n = 0, nterms 
            do m = 0, n 
              cquadr1 = (0d0, 0d0)
              cquadr2 = (0d0, 0d0)
              do i = 1, ntheta
                ctmp1 = (0d0, 0d0)
                ctmp2 = (0d0, 0d0)
                do j = 1, nphi
                  ctmp1 = ctmp1 + sphfgrid(j,i,n,m) * ws1grid(j,i)
                  ctmp2 = ctmp2 + sphfgrid(j,i,n,m) * ws2grid(j,i)
                enddo 
                cquadr1 = cquadr1 + ctmp1 * whts(i)
                cquadr2 = cquadr2 + ctmp2 * whts(i)
              enddo 
              cdrop1_chrg_LHS(n,m) = cquadr1 * 2 * pi / nphi
              cdrop2_chrg_LHS(n,m) = cquadr2 * 2 * pi / nphi
              if (m > 0) then
                cdrop1_chrg_LHS(n,-m) = dconjg(cdrop1_chrg_LHS(n,m))
                cdrop2_chrg_LHS(n,-m) = dconjg(cdrop2_chrg_LHS(n,m))
              endif
            enddo 
          enddo  
  !ccccccccccccccccc build linear system ccccccccccccccccccc
  allocate(cLHS_mat(N2,N2))
  allocate(cRHS_col(N2))

  !!!! the LHS matrix made up of 9 blocks
  !     N21*N21 + N21*N21 + N21*2
  !     N21*N21 + N21*N21 + N21*2
  !     2*N21  + 2*N21  + 2*2  
  !!!!

  !!!!!! net charge eqn  !!!!
          k = 1
          do n = 0, nterms 
            do m = -n, n 
              cLHS_mat(N2-1, k) = cdrop1_chrg_LHS(n,m) 
              cLHS_mat(N2, k + N21) = cdrop2_chrg_LHS(n,m) 
              cLHS_mat(N2-1, k + N21) = (0d0,0d0)
              cLHS_mat(N2, k) = (0d0,0d0)
              k = k + 1
            enddo 
          enddo  
          !!! zero net charge on drops
          cRHS_col(N2 - 1) = (0d0,0d0)
          cRHS_col(N2) = (0d0,0d0)

       j = 1
       do n = 0,nterms !! start 0
         do m = -n, n    !!! real, change to 0
           i = 1
           do n1 = 0,nterms !! start 0
             do m1 = -n1, n1   !!! real, change to 0
               cLHS_mat(i, j)=cls1d1mpoles(n1,m1,n,m)
               cLHS_mat(i, j + N21)=cls1d2mpoles(n1,m1,n,m)
               cLHS_mat(i + N21, j)=cls2d1mpoles(n1,m1,n,m)
               cLHS_mat(i + N21, j + N21)=cls2d2mpoles(n1,m1,n,m) 
  !        if (i == 1) then
        ! write(9401,*) i, j, n1,m1,n,m, cls1d1mpoles(n1,m1,n,m)!, 
  !     $ cls1d2mpoles(n1,m1,n,m),cls2d1mpoles(n1,m1,n,m),
  !     $ cls2d2mpoles(n1,m1,n,m)
  !        endif
               i = i + 1    
             enddo
           enddo
           j = j + 1
         enddo
       enddo    
       
       !! initialize the rightmost 2 cols of LHS
       do i = 1, N2
         cLHS_mat(i,N2-1) = (0d0,0d0)
         cLHS_mat(i,N2) = (0d0,0d0)
       enddo
       !! terms for the potential of conducting drops
       cLHS_mat(1,N2-1) = (1d0,0d0)
       cLHS_mat(N21+1,N2) = (1d0,0d0)


       !write(14,*) 'Y proj of E_inf*ngrid, ALL B terms in A * X = B'

       i = 1
       do n = 0, nterms !! start 0
           do m = -n, n   !!! real, change to 0           
               cRHS_col(i) = Pinf1mpole(n, m)
               cRHS_col(i + N21) = Pinf2mpole(n, m)
              ! write(14,*) 'n',n,'m',m,'i=',i,cap_Einf(i)
               i = i + 1
           enddo
       enddo

!!! screen tiny terms       
       tol_squ_mat = 1d-32

       write(21,*) "istep=",istep,"only nonzero B terms, in A * X = B","N21=",N21
       do i = 1, N2
!         if (abs(cRHS_col(i)*dconjg(cRHS_col(i))).lt.tol_squ_mat) then
!           cRHS_col(i) = 0d0
!         else
           write(21,*)  'i=',i, cRHS_col(i)
!         endif
       enddo
  
        !write(19,*)  " all A terms, in A * X = B"   
    !  write(20,*)  "istep=",istep, " only nonzero A terms, in A * X = B"  
      write(201,*)  "istep=",istep, " only nonzero imaginary of A terms, in A * X = B"        
        do k = 1, N2
          do l = 1, N2
            !write(19,*) 'i=',k,'j=',l, cap_lpl(k,l)
!            tmp = abs(cLHS_mat(k,l)*dconjg(cLHS_mat(k,l)))
!            if (tmp.lt.tol_squ_mat) then 
!              cLHS_mat(k,l) = 0d0
!            else
!           !  write(20,*) 'i=',k,'j=',l, cLHS_mat(k,l)
!            endif
            if (dimag(cLHS_mat(k,l)).gt.1d-13) then
             write(201,*) 'i=',k,'j=',l, cLHS_mat(k,l)
            endif
          enddo
        enddo
!!!!!!!!  

       !!!! special case for the first step sphere shape
       !! avoid all zero entry on the first line of linear system (A and B) 
  !     if (istep .eq. 1) then 
  !
  !       do j = 1, N2
  !        cLHS_mat(1,j) = (0d0, 0d0)
  !        cLHS_mat(1+N2/2,j) = (0d0, 0d0)
  !       enddo
  !       cLHS_mat(1,1) = (1d0, 0d0)
  !       cLHS_mat(1+N2/2,1+N2/2) = (1d0, 0d0)
  !     endif                 
        
  !cccccccccccc
       allocate(ipiv(N2))         !! n start with 0            
       CALL zgesv(N2,NRHS,cLHS_mat,N2,IPIV,cRHS_col,N2,INFO)      
  !      CALL zgesv(N2-1,NRHS,cap_lpl,N2-1,
  !     $ IPIV,cap_Einf,N2-1,INFO)      
       write(116,*) "istep=",istep,'linear solver zgesv solution:','INFO=',INFO
       write(216,*) "istep=",istep,'linear solver zgesv solution:','INFO=',INFO
       i = 1
       do n = 0, nterms !! start 0
           do m = -n, n  !!! real, change to 0   
  
               if (m .ge. 0) then         
               
               en1mpole(n, m) = cRHS_col(i)   
               en2mpole(n, m) = cRHS_col(i + N21)      
               
               write(116,*) 'i=',i,'n',n,'m',m,cRHS_col(i)
               write(216,*) 'i=',i,'n',n,'m',m,cRHS_col(i+ N21)
               tmp = abs(en1mpole(n, m)*dconjg(en1mpole(n, m)))
               if (tmp.lt.tol_squ_mat) then
                 en1mpole(n, m) = 0d0
               endif
               tmp = abs(en2mpole(n, m)*dconjg(en2mpole(n, m)))
               if (tmp.lt.tol_squ_mat) then
                 en2mpole(n, m) = 0d0
               endif
  
               endif
               i = i + 1
           enddo
       enddo
       
       !! potentials obtained
       write(116,*) 'i=',N2-1,cRHS_col(N2-1)
       write(216,*) 'i=',N2,cRHS_col(N2)
  !      write(111,*) "n, m, en1mpole(n,m), en2mpole(n,m)"
  !      do n = 1, nterms
  !        do m = 0, n
  !          write(111,*) n, m, en1mpole(n,m), en2mpole(n,m)
  !       enddo
  !      enddo
  !ccccc spherical harmonic transform ccccccc
       allocate(En1(nphi, ntheta),En2(nphi, ntheta))
       write(17,*) "istep=",istep,"'j=',j,'i=',i, En1(j,i),En2(j,i)"
       do i = 1,ntheta
           do j = 1, nphi
             ctmp1 = (0d0,0d0)
             ctmp2 = (0d0,0d0)
             do n = 0, nterms   !! start 0
               m = 0
               ctmp1 = ctmp1 + en1mpole(n, m)*sphfgrid(j,i,n,m) 
               ctmp2 = ctmp2 + en2mpole(n, m)*sphfgrid(j,i,n,m)                    
               do m = 1, n    !!! real, change to 0
                   ctmp1 = ctmp1 + 2*en1mpole(n, m)*sphfgrid(j,i,n,m) 
                   ctmp2 = ctmp2 + 2*en2mpole(n, m)*sphfgrid(j,i,n,m)                
               enddo
             enddo 
             En1(j,i) = dble(ctmp1)
             En2(j,i) = dble(ctmp2)
             write(17,*) 'j=',j,'i=',i, En1(j,i),En2(j,i)
           enddo
       enddo
  !cccccccccc total force mpole cccccccccc
  !c          ( En^2/2*Ca + 2*H ) * ngrid
  !ccccccccccccccccccccccccccccccccccccccc    
        write(*,*) "start: evaluate velocity"
  
       allocate(fx1(nphi,ntheta), fy1(nphi,ntheta), fz1(nphi,ntheta), &
        fx2(nphi,ntheta), fy2(nphi,ntheta), fz2(nphi,ntheta))
  
       Caflow = 0.0d0
       write(121,*) "istep=",istep,"i,j,tmp1,fx1(j,i),fy1(j,i),fz1(j,i)"
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
           write(121,*) i,j,tmp1,fx1(j,i),fy1(j,i),fz1(j,i)
         enddo
       enddo
  
       allocate(fx1mpole(0:nterms,0:nterms), fy1mpole(0:nterms,0:nterms),&
        fz1mpole(0:nterms,0:nterms),fx2mpole(0:nterms,0:nterms),&
        fy2mpole(0:nterms,0:nterms),fz2mpole(0:nterms,0:nterms))
       
       call sphtrans_fwd_real_rsv(nterms,fx1mpole,nphi,ntheta,fx1,&
            ctheta,whts,ynms,dwsave)
       call sphtrans_fwd_real_rsv(nterms,fy1mpole,nphi,ntheta,fy1,&
            ctheta,whts,ynms,dwsave)    
       call sphtrans_fwd_real_rsv(nterms,fz1mpole,nphi,ntheta,fz1,&
            ctheta,whts,ynms,dwsave)           
       call sphtrans_fwd_real_rsv(nterms,fx2mpole,nphi,ntheta,fx2,&
            ctheta,whts,ynms,dwsave)
       call sphtrans_fwd_real_rsv(nterms,fy2mpole,nphi,ntheta,fy2,&
            ctheta,whts,ynms,dwsave)    
       call sphtrans_fwd_real_rsv(nterms,fz2mpole,nphi,ntheta,fz2,&
            ctheta,whts,ynms,dwsave)           
  !cccccccccccccccccccccccccccccccccccccccc
  
  !cccccc rotate force components mpole ccccccccccc
  
       allocate(rfx1grids(nphi,ntheta,nrot,nbeta), &
        rfy1grids(nphi,ntheta,nrot,nbeta), &
        rfz1grids(nphi,ntheta,nrot,nbeta), &
        rfx2grids(nphi,ntheta,nrot,nbeta), &
        rfy2grids(nphi,ntheta,nrot,nbeta), &
        rfz2grids(nphi,ntheta,nrot,nbeta))
  
       call rotgrid_fsr_real(nterms,fx1mpole,nphi,ntheta, &
        nbeta,theta,nrot,rfx1grids, &
        rotmat,ctheta,ynms,dwsave)
       call rotgrid_fsr_real(nterms,fy1mpole,nphi,ntheta, &
        nbeta,theta,nrot,rfy1grids, &
        rotmat,ctheta,ynms,dwsave)
       call rotgrid_fsr_real(nterms,fz1mpole,nphi,ntheta, &
        nbeta,theta,nrot,rfz1grids, &
        rotmat,ctheta,ynms,dwsave)            
       call rotgrid_fsr_real(nterms,fx2mpole,nphi,ntheta, &
        nbeta,theta,nrot,rfx2grids, &
        rotmat,ctheta,ynms,dwsave)
       call rotgrid_fsr_real(nterms,fy2mpole,nphi,ntheta, &
        nbeta,theta,nrot,rfy2grids, &
        rotmat,ctheta,ynms,dwsave)
       call rotgrid_fsr_real(nterms,fz2mpole,nphi,ntheta, &
        nbeta,theta,nrot,rfz2grids, &
        rotmat,ctheta,ynms,dwsave)    
  !cccccccccccccccccccccccccccccccccccccccccccccc
       
       allocate(us1d1grid(nphi,ntheta,3),us2d2grid(nphi,ntheta,3), &
        us1d2grid(nphi,ntheta,3),us2d1grid(nphi,ntheta,3))
  !cc    special gauss-legendre quadrature(lapalce) for all rotated grids
           !write(34,*) 'irot ','jbeta ', 'special g-l quadr stk'
           do kk = 1,nrot
             do k = 1,nbeta
               call integral_eval_ufield(ntheta,nphi, &
                rx1grids(:,:,kk,k), ry1grids(:,:,kk,k), rz1grids(:,:,kk,k),  &
                rws1grids(:,:,kk,k), rfx1grids(:,:,kk,k), rfy1grids(:,:,kk,k),  &
                rfz1grids(:,:,kk,k), swhts,  &
                x1grid(kk,k), y1grid(kk,k), z1grid(kk,k),  &
                quadr1, quadr2, quadr3)
               us1d1grid(kk,k,1) = quadr1 
               us1d1grid(kk,k,2) = quadr2
               us1d1grid(kk,k,3) = quadr3 
               !write(34,*) kk, k, us1d1grid(kk,k,1), us1d1grid(kk,k,2),
    !  $  us1d1grid(kk,k,3)
             enddo
           enddo 
  
           do kk = 1,nrot
             do k = 1,nbeta
               call integral_eval_ufield(ntheta,nphi, &
                rx2grids(:,:,kk,k), ry2grids(:,:,kk,k), rz2grids(:,:,kk,k),  &
                rws2grids(:,:,kk,k), rfx2grids(:,:,kk,k), rfy2grids(:,:,kk,k),  &
                rfz2grids(:,:,kk,k), swhts,  &
                x2grid(kk,k), y2grid(kk,k), z2grid(kk,k),  &
                quadr1, quadr2, quadr3)
               us2d2grid(kk,k,1) = quadr1
               us2d2grid(kk,k,2) = quadr2 
               us2d2grid(kk,k,3) = quadr3
              ! write(34,*) kk, k, us2d2grid(kk,k,1), us2d2grid(kk,k,2),
     ! $  us2d2grid(kk,k,3)
             enddo
           enddo 
  
       allocate(u1mpoles(0:nterms,0:nterms,3), &
        u2mpoles(0:nterms,0:nterms,3))
  
       do i = 1, 3
              call sphtrans_fwd_real_rsv(nterms,u1mpoles(:,:,i), &
                  nphi,ntheta,us1d1grid(:,:,i),ctheta,whts,ynms,dwsave)
              call sphtrans_fwd_real_rsv(nterms,u2mpoles(:,:,i), &
                  nphi,ntheta,us2d2grid(:,:,i),ctheta,whts,ynms,dwsave)      
       enddo
  
  
       allocate(fx1igrid(iphi,itheta), fy1igrid(iphi,itheta),&
        fz1igrid(iphi,itheta), fx2igrid(iphi,itheta), &
        fy2igrid(iphi,itheta), fz2igrid(iphi,itheta))
  
       allocate(En1igrid(iphi,itheta), En2igrid(iphi,itheta))
  
       call nearSing_upsamp_En(iup, iterms, itheta, iphi, &
         en1mpole, en2mpole, En1igrid, En2igrid)
       
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
               call close_dist_flag(iflag_cd, x1grid(kk,k),y1grid(kk,k),z1grid(kk,k), &
               x2grid, y2grid, z2grid, x2mpole, y2mpole, z2mpole, &
               dist_cl, sphf_s, xs, ys, zs) !! x,y,zs are the closest point on surface
  
               if (dist_cl .ge. h_far) then 
  
                 call integral_eval_ufield(ntheta, nphi, &
                           x2grid, y2grid, z2grid,  &
                           ws2grid, fx2, fy2, fz2, whts, &    !! regular whts 
                           x1grid(kk,k), y1grid(kk,k), z1grid(kk,k),  &
                           quadr1, quadr2, quadr3)
  
               else 
                 !! call near_singlar
                 call nearSing_integral_ufield(iterms, itheta, iphi, &
                           x2igrid, y2igrid, z2igrid, ws2igrid,  &
                           fx2igrid, fy2igrid, fz2igrid, whts_i, &
                           u2mpoles, dist_cl, sphf_s, & !! interpolate singluar integral
                           xs, ys, zs, x1grid(kk,k),y1grid(kk,k),z1grid(kk,k), &
                           quadr1, quadr2, quadr3)
  
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
                  
                 call close_dist_flag(iflag_cd, x2grid(kk,k),y2grid(kk,k),z2grid(kk,k), &
                           x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole, &
                           dist_cl, sphf_s, xs, ys, zs) !! x,y,zs are the closest point on surface
   
               if (dist_cl .ge. h_far) then 
  
                 call integral_eval_ufield(ntheta, nphi, &
                           x1grid, y1grid, z1grid,  &
                           ws1grid, fx1, fy1, fz1, whts, &    !! regular whts 
                           x2grid(kk,k), y2grid(kk,k), z2grid(kk,k),  &
                           quadr1, quadr2, quadr3)
  
               else 
                 !! call near_singlar
                 call nearSing_integral_ufield(iterms, itheta, iphi,  &
                           x1igrid, y1igrid, z1igrid, ws1igrid,  &
                           fx1igrid, fy1igrid, fz1igrid, whts_i, &
                           u1mpoles, dist_cl, sphf_s,  & !! interpolate singluar integral 
                           xs, ys, zs, x2grid(kk,k), y2grid(kk,k), z2grid(kk,k), &
                           quadr1, quadr2, quadr3)
  
              endif
              us2d1grid(kk,k,1) = quadr1
              us2d1grid(kk,k,2) = quadr2 
              us2d1grid(kk,k,3) = quadr3
             enddo
           enddo    
           write(135,*) "istep=",istep,"i,j,k, us1d1grid(j,i,k), us1d2grid(j,i,k)"
           write(235,*) "istep=",istep,"i,j,k, us2d1grid(j,i,k), us2d2grid(j,i,k)"
           do i = 1, ntheta
             do j = 1, nphi
               do k = 1,3
              write(135,*) i,j,k, us1d1grid(j,i,k), us1d2grid(j,i,k)
              write(235,*) i,j,k, us2d1grid(j,i,k), us2d2grid(j,i,k)
               u1grid(j,i,k) = us1d1grid(j,i,k) + us1d2grid(j,i,k)
               u2grid(j,i,k) = us2d1grid(j,i,k) + us2d2grid(j,i,k)
               enddo
             enddo
           enddo
  
         write(*,*) "eval velocity done"
  
  !!!<       electrical total force on drop 1
         ele_tot_fx = 0d0
         ele_tot_fy = 0d0
         ele_tot_fz = 0d0
         charge_net = 0d0
         do i = 1, ntheta
           tmpx = 0d0
           tmpy = 0d0
           tmpz = 0d0
           tmpc = 0d0
           do j = 1, nphi            
             tmp = En1(j,i)**2/2d0*Ca
             tmpx = tmpx + tmp*n1grid(j,i,1)*ws1grid(j,i)
             tmpy = tmpy + tmp*n1grid(j,i,2)*ws1grid(j,i)
             tmpz = tmpz + tmp*n1grid(j,i,3)*ws1grid(j,i)
             tmpc = tmpc + En1(j,i)*ws1grid(j,i) 
           enddo
           ele_tot_fx = ele_tot_fx + tmpx * whts(i)
           ele_tot_fy = ele_tot_fy + tmpy * whts(i)
           ele_tot_fz = ele_tot_fz + tmpz * whts(i)
           charge_net = charge_net + tmpc * whts(i)
         enddo
         ele_tot_fx = ele_tot_fx/nphi*(2*pi)
         ele_tot_fy = ele_tot_fy/nphi*(2*pi)
         ele_tot_fz = ele_tot_fz/nphi*(2*pi)
         ele_tot_fmag = ele_tot_fx**2 + ele_tot_fy**2 + ele_tot_fz**2
         ele_tot_fmag = dsqrt(ele_tot_fmag)
         charge_net = charge_net/nphi*(2*pi)
         write(109,*) ele_tot_fx, ele_tot_fy, ele_tot_fz, ele_tot_fmag, charge_net
  
         tot_fx = 0d0
         tot_fy = 0d0
         tot_fz = 0d0
         do i = 1, ntheta
           tmpx = 0d0
           tmpy = 0d0
           tmpz = 0d0
           do j = 1, nphi            
             tmp = 2d0*h1grid(j,i)
             tmpx = tmpx + tmp*n1grid(j,i,1)*ws1grid(j,i)
             tmpy = tmpy + tmp*n1grid(j,i,2)*ws1grid(j,i)
             tmpz = tmpz + tmp*n1grid(j,i,3)*ws1grid(j,i)
           enddo
           tot_fx = tot_fx + tmpx * whts(i)
           tot_fy = tot_fy + tmpy * whts(i)
           tot_fz = tot_fz + tmpz * whts(i)
         enddo
         tot_fx = tot_fx/nphi*(2*pi)
         tot_fy = tot_fy/nphi*(2*pi)
         tot_fz = tot_fz/nphi*(2*pi)
         tot_fmag = tot_fx**2 + tot_fy**2 + tot_fz**2
         tot_fmag = dsqrt(tot_fmag)
         write(108,*) tot_fx, tot_fy, tot_fz, tot_fmag  
  !!!>

       deallocate(x1igrid, y1igrid, z1igrid, ws1igrid, h1igrid, &
        x2igrid, y2igrid, z2igrid, ws2igrid, h2igrid, whts_i)
       deallocate(n1igrid, n2igrid, cpinf1grid,cpinf2grid)
       deallocate(rx1grids, ry1grids, &
        rz1grids, rws1grids, rx2grids, ry2grids, rz2grids, rws2grids)
      !!deallocate(f1mpoles,f2mpoles)
       deallocate(cls1d1grids, cls1d2grids, cls2d1grids, cls2d2grids)
       deallocate(cls1d1mpoles,cls1d2mpoles, cls2d1mpoles, cls2d2mpoles)
       deallocate(Pinf1mpole, Pinf2mpole, cLHS_mat, cRHS_col, ipiv)
       deallocate(En1, En2, fx1, fy1, fz1, fx2, fy2, fz2)
       deallocate(fx1mpole, fy1mpole, fz1mpole, fx2mpole, fy2mpole,  &
        fz2mpole)
       deallocate(rfx1grids, rfy1grids, rfz1grids, rfx2grids, rfy2grids,  &
        rfz2grids)
       deallocate(us1d1grid, us2d1grid, us1d2grid, us2d2grid)
       deallocate(u1mpoles, u2mpoles)
       deallocate(fx1igrid, fy1igrid, fz1igrid, fx2igrid, fy2igrid,  &
        fz2igrid, En1igrid, En2igrid)
  
       end SUBROUTINE
  
  
  !c  Our definition of complex spherical harmonics is
  !c
  !c  Ynm(theta,phi)= sqrt(2n+1) sqrt((n-m)!/(n+m)!) 
  !c                  Pnm(cos theta) e^(im phi), 
  !c  Yn,-m(theta,phi) = sqrt(2n+1) sqrt((n-m)!/(n+m)!) 
  !c                  Pnm(cos theta) e^(-im phi),   for m >= 0.
  !c       
  !c  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.
  
  
  
  
