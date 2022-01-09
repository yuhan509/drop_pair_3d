
c
c
c     
      subroutine shf_curv_up(iterms, itheta, iphi, iup, 
     $ ximpole, yimpole, zimpole,
     $ higrid, wigrid, nigrid)
     
      use mod_upsamp

      implicit real *8 (a-h,o-z)
      integer :: iterms, itheta, iphi, iup
      complex *16 :: ximpole(0:iterms,0:iterms),
     $ yimpole(0:iterms,0:iterms),zimpole(0:iterms,0:iterms)
      real *8 :: higrid(iphi, itheta), wigrid(iphi, itheta),
     $ nigrid(iphi, itheta,3)

      complex *16, ALLOCATABLE, dimension(:,:,:,:) ::  !! change to real *8
     $ dydph, dydth, d2ydth2, d2ydph2, d2ydthdph  
      real* 8, ALLOCATABLE, dimension(:,:,:):: dxdth,d2xdth,
     $ dxdph,d2xdph,d2xdthdph
      real* 8, ALLOCATABLE:: Geom_E(:,:),Geom_F(:,:),Geom_G(:,:)
      real* 8, ALLOCATABLE:: Geom_L(:,:),Geom_M(:,:),Geom_N(:,:)   

      allocate(dxdth(iphi,itheta,3))
      allocate(dxdph(iphi,itheta,3))
      allocate(d2xdth(iphi,itheta,3))
      allocate(d2xdph(iphi,itheta,3))
      allocate(d2xdthdph(iphi,itheta,3))
      allocate(Geom_E(iphi,itheta))
      allocate(Geom_F(iphi,itheta))
      allocate(Geom_G(iphi,itheta))
      allocate(Geom_L(iphi,itheta))
      allocate(Geom_M(iphi,itheta))
      allocate(Geom_N(iphi,itheta))

      allocate(dydth(iphi,itheta,0:iterms,-iterms:iterms))
      allocate(dydph(iphi,itheta,0:iterms,-iterms:iterms))
      allocate(d2ydth2(iphi,itheta,0:iterms,-iterms:iterms))
      allocate(d2ydthdph(iphi,itheta,0:iterms,-iterms:iterms))
      allocate(d2ydph2(iphi,itheta,0:iterms,-iterms:iterms))

          do n = 0, iterms   !! iterms to nterms
            do m = 0, n   !!! real, change to zero
              do i = 1, itheta
                do j = 1, iphi
                  dydph(j,i,n,m) = up_dydph(j,i,n,m,iup)
                  dydth(j,i,n,m) = up_dydth(j,i,n,m,iup) 
                  d2ydph2(j,i,n,m) = up_d2ydph2(j,i,n,m,iup)
                  d2ydth2(j,i,n,m) = up_d2ydth2(j,i,n,m,iup)
                  d2ydthdph(j,i,n,m) = up_d2ydthdph(j,i,n,m,iup)                 
                enddo
              enddo
            enddo
          enddo 
      
      do n = 0, iterms
        do m = 0, n
          if (imag(ximpole(n,m)) .gt. 1d-14) then
            write(113,*) "x",n,m,ximpole(n,m)
          endif
          if (imag(yimpole(n,m)) .gt. 1d-14) then
            write(113,*) "y",n,m,yimpole(n,m)
          endif
          if (imag(zimpole(n,m)) .gt. 1d-14) then
            write(113,*) "z",n,m,zimpole(n,m)
         endif
        enddo
      enddo
      do i=1,itheta
        do j=1,iphi
          do k = 1,3
cc        first order derivatives              
          dxdth(j,i,k) = 0d0                
          dxdph(j,i,k) = 0d0 
cc        second order derivatives    
          d2xdth(j,i,k) = 0d0
          d2xdph(j,i,k) = 0d0
          d2xdthdph(j,i,k) = 0d0            
          enddo
   
          do l = 0,iterms
            m = 0
cc        first order derivatives
        dxdth(j,i,1)=dxdth(j,i,1)+dble(ximpole(l,m)*dydth(j,i,l,m))
        dxdth(j,i,2)=dxdth(j,i,2)+dble(yimpole(l,m)*dydth(j,i,l,m))
        dxdth(j,i,3)=dxdth(j,i,3)+dble(zimpole(l,m)*dydth(j,i,l,m))
          
        dxdph(j,i,1)=dxdph(j,i,1)+dble(ximpole(l,m)*dydph(j,i,l,m))
        dxdph(j,i,2)=dxdph(j,i,2)+dble(yimpole(l,m)*dydph(j,i,l,m))

        dxdph(j,i,3)=dxdph(j,i,3)+dble(zimpole(l,m)*dydph(j,i,l,m))
cc        second order derivatives        
            d2xdth(j,i,1)=d2xdth(j,i,1)+
     $ dble(ximpole(l,m)*d2ydth2(j,i,l,m))
            d2xdth(j,i,2)=d2xdth(j,i,2)+
     $ dble(yimpole(l,m)*d2ydth2(j,i,l,m))
            d2xdth(j,i,3)=d2xdth(j,i,3)+
     $ dble(zimpole(l,m)*d2ydth2(j,i,l,m))
     
            d2xdph(j,i,1)=d2xdph(j,i,1)+
     $ dble(ximpole(l,m)*d2ydph2(j,i,l,m))
            d2xdph(j,i,2)=d2xdph(j,i,2)+
     $ dble(yimpole(l,m)*d2ydph2(j,i,l,m))
            d2xdph(j,i,3)=d2xdph(j,i,3)+
     $ dble(zimpole(l,m)*d2ydph2(j,i,l,m))

            d2xdthdph(j,i,1)=d2xdthdph(j,i,1)+            
     $ dble(ximpole(l,m)*d2ydthdph(j,i,l,m))
            d2xdthdph(j,i,2)=d2xdthdph(j,i,2)+
     $ dble(yimpole(l,m)*d2ydthdph(j,i,l,m))
            d2xdthdph(j,i,3)=d2xdthdph(j,i,3)+
     $ dble(zimpole(l,m)*d2ydthdph(j,i,l,m))    
!      write(42,*) j,i,l,m,"d2ydth2",d2ydth2(j,i,l,m)
!      write(42,*) j,i,l,m,"d2ydph2", d2ydph2(j,i,l,m)
!      write(42,*) j,i,l,m,"d2ydthdph",d2ydthdph(j,i,l,m)           
            do m = 1, l 
cc        first order derivatives
        dxdth(j,i,1)=dxdth(j,i,1)+ 2*dble(ximpole(l,m)*dydth(j,i,l,m))
        dxdth(j,i,2)=dxdth(j,i,2)+ 2*dble(yimpole(l,m)*dydth(j,i,l,m))
        dxdth(j,i,3)=dxdth(j,i,3)+ 2*dble(zimpole(l,m)*dydth(j,i,l,m))

        dxdph(j,i,1)=dxdph(j,i,1)+ 2*dble(ximpole(l,m)*dydph(j,i,l,m))
        dxdph(j,i,2)=dxdph(j,i,2)+ 2*dble(yimpole(l,m)*dydph(j,i,l,m))
        dxdph(j,i,3)=dxdph(j,i,3)+ 2*dble(zimpole(l,m)*dydph(j,i,l,m))
cc        second order derivatives        
            d2xdth(j,i,1)=d2xdth(j,i,1)+
     $ 2 * dble(ximpole(l,m)*d2ydth2(j,i,l,m))
            d2xdth(j,i,2)=d2xdth(j,i,2)+
     $ 2 * dble(yimpole(l,m)*d2ydth2(j,i,l,m))
            d2xdth(j,i,3)=d2xdth(j,i,3)+
     $ 2 * dble(zimpole(l,m)*d2ydth2(j,i,l,m))
     
            d2xdph(j,i,1)=d2xdph(j,i,1)+
     $ 2 * dble(ximpole(l,m)*d2ydph2(j,i,l,m))
            d2xdph(j,i,2)=d2xdph(j,i,2)+
     $ 2 * dble(yimpole(l,m)*d2ydph2(j,i,l,m))
            d2xdph(j,i,3)=d2xdph(j,i,3)+
     $ 2 * dble(zimpole(l,m)*d2ydph2(j,i,l,m))

            d2xdthdph(j,i,1)=d2xdthdph(j,i,1)+            
     $ 2 * dble(ximpole(l,m)*d2ydthdph(j,i,l,m))
            d2xdthdph(j,i,2)=d2xdthdph(j,i,2)+
     $ 2 * dble(yimpole(l,m)*d2ydthdph(j,i,l,m))
            d2xdthdph(j,i,3)=d2xdthdph(j,i,3)+
     $ 2 * dble(zimpole(l,m)*d2ydthdph(j,i,l,m))    
!          write(42,*) j,i,l,m,"d2ydth2",d2ydth2(j,i,l,m)
!          write(42,*) j,i,l,m,"d2ydph2", d2ydph2(j,i,l,m)
!          write(42,*) j,i,l,m,"d2ydthdph",d2ydthdph(j,i,l,m)     
            enddo
          enddo

!         write(41,*) j,i,"dxdth",dxdth(j,i,1),dxdth(j,i,2),dxdth(j,i,3)
!         write(41,*) j,i,"dxdph",dxdph(j,i,1),dxdph(j,i,2),dxdph(j,i,3)
!      write(41,*) j,i,"d2xdth",d2xdth(j,i,1),d2xdth(j,i,2),d2xdth(j,i,3)
!      write(41,*) j,i,"d2xdph",d2xdph(j,i,1),d2xdph(j,i,2),d2xdph(j,i,3)
!      write(41,*) j,i,"d2xdthdph",
!    $ d2xdthdph(j,i,1),d2xdthdph(j,i,2),d2xdthdph(j,i,3)

          Geom_E(j,i)=dxdth(j,i,1)*dxdth(j,i,1)
     $    + dxdth(j,i,2)*dxdth(j,i,2) + dxdth(j,i,3)*dxdth(j,i,3)
          Geom_F(j,i)=dxdth(j,i,1)*dxdph(j,i,1)
     $    + dxdth(j,i,2)*dxdph(j,i,2) + dxdth(j,i,3)*dxdph(j,i,3)
          Geom_G(j,i)=dxdph(j,i,1)*dxdph(j,i,1)
     $    + dxdph(j,i,2)*dxdph(j,i,2) + dxdph(j,i,3)*dxdph(j,i,3) 
          wigrid(j,i)=Geom_E(j,i)*Geom_G(j,i) - Geom_F(j,i)**2
          wigrid(j,i)=dsqrt(wigrid(j,i))

      nigrid(j,i,1)=dxdth(j,i,2)*dxdph(j,i,3)-dxdth(j,i,3)*dxdph(j,i,2)
      nigrid(j,i,2)=dxdth(j,i,3)*dxdph(j,i,1)-dxdth(j,i,1)*dxdph(j,i,3)
      nigrid(j,i,3)=dxdth(j,i,1)*dxdph(j,i,2)-dxdth(j,i,2)*dxdph(j,i,1)
        nigrid(j,i,1) =nigrid(j,i,1) / wigrid(j,i)
        nigrid(j,i,2) =nigrid(j,i,2) / wigrid(j,i)
        nigrid(j,i,3) =nigrid(j,i,3) / wigrid(j,i)

          Geom_L(j,i) = d2xdth(j,i,1) *nigrid(j,i,1) 
     $ + d2xdth(j,i,2) *nigrid(j,i,2) + d2xdth(j,i,3) *nigrid(j,i,3)
          Geom_M(j,i) = d2xdthdph(j,i,1) * nigrid(j,i,1) 
     $ + d2xdthdph(j,i,2)*nigrid(j,i,2)+d2xdthdph(j,i,3)*nigrid(j,i,3)
          Geom_N(j,i) = d2xdph(j,i,1) * nigrid(j,i,1) 
     $ + d2xdph(j,i,2) *nigrid(j,i,2) + d2xdph(j,i,3) *nigrid(j,i,3)
          higrid(j,i) = (Geom_E(j,i)*Geom_N(j,i)+Geom_G(j,i)*Geom_L(j,i)
     $ - 2d0*Geom_F(j,i)*Geom_M(j,i))/(2d0*wigrid(j,i)**2)
!          Geom_K(j,i) = (Geom_L(j,i)*Geom_N(j,i) - Geom_M(j,i)**2)
!     $ / (wigrid(j,i)**2)
        enddo
      enddo 

      END SUBROUTINE
c
c
c
      subroutine eval_curv_mpole_uphonly(xmpole, ympole, zmpole, 
     $ hgrid, wgrid, ngrid, hmpole)
      use mod_upsamp
      implicit real *8 (a-h,o-z)
      !! input
      complex *16 :: xmpole(0:nterms,0:nterms),
     $ ympole(0:nterms,0:nterms),zmpole(0:nterms,0:nterms)
      
      !! output
      real *8 :: hgrid(nphi, ntheta), wsgrid(nphi, ntheta), 
     $ ngrid(nphi, ntheta,3)
      complex *16 :: hmpole(0:nterms,0:nterms)
 
      integer :: iterms, itheta, iphi, iup
      complex *16, ALLOCATABLE :: himpole(:,:)
      real *8, ALLOCATABLE :: higrid(:,:),wigrid(:,:), nigrid(:,:,:)
      real* 8, ALLOCATABLE:: ctheta(:), whts(:), theta(:), ynms(:,:,:) 
      complex *16, ALLOCATABLE:: dwsave(:)    
      

      iup = 1
      write(*,*) "start of iup = 1"

      call shf_curv_up(nterms, ntheta, nphi, iup,
     $ xmpole, ympole, zmpole, 
     $ hgrid, wsgrid, ngrid)

      allocate(ctheta(ntheta), whts(ntheta), theta(ntheta))
      allocate(ynms(0:nterms,0:nterms,ntheta/2+1))
      allocate(dwsave(4*nphi+15))

      do i = 1, ntheta
        ctheta(i) = up_ctheta(i,iup) 
        whts(i) = up_whts(i,iup)
        theta(i) = up_theta(i,iup)
      enddo
      do j = 1, nphi*4+15
        dwsave(j) = up_dwsave(j,iup)
      enddo
      do n = 0, nterms
        do m = 0, n
          do i = 1, ntheta/2+1
             ynms(n,m,i) = up_ynms(n,m,i,iup)
          enddo
        enddo
      enddo 
      
      !! use wgrid/dsin(theta(i)) instead of wgrid

       do i = 1, ntheta
        do j = 1, nphi
         wsgrid(j,i) = wsgrid(j,i) / dsin(theta(i))
        enddo
       enddo

      call sphtrans_fwd_real(nterms,hmpole,nphi,ntheta,hgrid,
     $ ctheta,whts,ynms,dwsave)     

      write(*,*) "end of iup = 1"
      deallocate(ctheta, whts, ynms, dwsave)
      
      do iup = 2, maxup_crv 
       iterms = iup * nterms
       call grid_dim_init(iterms,itheta,iphi)

       allocate(himpole(0:iterms,0:iterms))
       allocate(higrid(iphi,itheta))
       allocate(wigrid(iphi,itheta))
       allocate(nigrid(iphi,itheta,3))
      
       write(*,*) "start evaluate higrid,wigrid,nigrid "

       call shf_curv_up(nterms, itheta, iphi, iup, 
     $ xmpole, ympole, zmpole, 
     $ higrid, wigrid, nigrid)
       
       write(*,*) "evaluate higrid,wigrid,nigrid done"
     
c       write(82,*) "iup =", iup
c       do i = 1,itheta
c        do j = 1, iphi
c          write(82,*) i,j,"H",higrid(j,i),"W",wigrid(j,i) 
c        enddo
c      enddo

       allocate(ctheta(itheta), whts(itheta))
       allocate(ynms(0:iterms,0:iterms,itheta/2+1))
       allocate(dwsave(4*iphi+15))
      
      !! load ctheta, whts, ynms from up_...
      do i = 1, itheta
        ctheta(i) = up_ctheta(i,iup) 
        whts(i) = up_whts(i,iup)
      enddo
      do j = 1, iphi*4+15
        dwsave(j) = up_dwsave(j,iup)
      enddo
      do n = 0, iterms
        do m = 0, n
          do i = 1, itheta/2+1
             ynms(n,m,i) = up_ynms(n,m,i,iup)
          enddo
        enddo
      enddo 
      
      call sphtrans_fwd_real_rsv(iterms,himpole,iphi,itheta,higrid,
     $ ctheta,whts,ynms,dwsave) 

      call eval_ups_tol_real(iterms, nterms, hmpole, himpole, res)

      !! trunc to order nterms hmpole
      do n = 0,nterms
        do m = 0,n
          hmpole(n, m) = himpole(n, m)    
        enddo
      enddo  

      !! if res < tol do the follow then break
     
      !! break the loop

      deallocate(ctheta, whts, ynms, dwsave)
      deallocate(higrid, wigrid, nigrid)
      deallocate(himpole)
      enddo

      allocate(ctheta(ntheta), whts(ntheta))
      allocate(ynms(0:nterms,0:nterms,ntheta/2+1))
      allocate(dwsave(4*nphi+15))

      iup = 1  !! dont forget reset iup
      do i = 1, ntheta 
        ctheta(i) = up_ctheta(i,iup) 
        whts(i) = up_whts(i,iup)
        theta(i) = up_theta(i,iup)
      enddo
      do j = 1, nphi*4+15
        dwsave(j) = up_dwsave(j,iup)
      enddo
      do n = 0, nterms
        do m = 0, n
          do i = 1, ntheta/2+1
             ynms(n,m,i) = up_ynms(n,m,i,iup)
          enddo
        enddo
      enddo 

      call sphtrans_real(nterms,hmpole,nphi,ntheta,hgrid,
     $ ctheta,ynms,dwsave)    
     
      END SUBROUTINE
c
c
c
c
      subroutine eval_curv_mpole_upAll(idrop, istep,
     $ xmpole, ympole, zmpole, 
     $  hgrid, wsgrid, ngrid, hmpole, wsmpole, nmpole)
      use mod_upsamp
      implicit real *8 (a-h,o-z)

      integer :: idrop
      complex *16 :: xmpole(0:nterms,0:nterms),
     $ ympole(0:nterms,0:nterms), zmpole(0:nterms,0:nterms)
      
      !! output
      real *8 :: hgrid(nphi, ntheta), wsgrid(nphi, ntheta), 
     $ ngrid(nphi, ntheta,3)
      complex *16 ::hmpole(0:nterms,0:nterms),
     $ wsmpole(0:nterms,0:nterms), nmpole(0:nterms,0:nterms,3)
 
      integer :: iterms, itheta, iphi, iup
      complex *16, ALLOCATABLE :: himpole(:,:), 
     $ wsimpole(:,:), nimpole(:,:,:) 
      real *8, ALLOCATABLE :: higrid(:,:),wsigrid(:,:), nigrid(:,:,:)
      real* 8, ALLOCATABLE:: ctheta(:), theta(:), whts(:), ynms(:,:,:) 
      complex *16, ALLOCATABLE:: dwsave(:)    
      
      iup = 1
      write(*,*) "start of iup = 1"

      call shf_curv_up(nterms,ntheta,nphi,iup,
     $ xmpole,ympole,zmpole, 
     $ hgrid,wsgrid,ngrid)

      allocate(ctheta(ntheta), whts(ntheta), theta(ntheta))
      allocate(ynms(0:nterms,0:nterms,ntheta/2+1))
      allocate(dwsave(4*nphi+15))

      do i = 1, ntheta
        ctheta(i) = up_ctheta(i,iup) 
        whts(i) = up_whts(i,iup)
        theta(i) = up_theta(i,iup)
      enddo
      do j = 1, nphi*4+15
        dwsave(j) = up_dwsave(j,iup)
      enddo
      do n = 0, nterms
        do m = 0, n
          do i = 1, ntheta/2+1
             ynms(n,m,i) = up_ynms(n,m,i,iup)
          enddo
        enddo
      enddo 
      
       do i = 1, ntheta
         do j = 1, nphi
          wsgrid(j,i) = wsgrid(j,i) / dsin(theta(i))
         enddo
       enddo

      call sphtrans_fwd_real_rsv(nterms,hmpole,nphi,ntheta,hgrid,
     $ ctheta,whts,ynms,dwsave)     

      
      deallocate(ctheta, theta, whts, ynms, dwsave)
      write(*,*) "end of iup = 1"

      do iup = 2, maxup_crv
       iterms = iup * nterms
       call grid_dim_init(iterms,itheta,iphi)
       
       allocate(himpole(0:iterms,0:iterms))
       allocate(wsimpole(0:iterms,0:iterms))
       allocate(nimpole(0:iterms,0:iterms,3))
       allocate(higrid(iphi,itheta))
       allocate(wsigrid(iphi,itheta))
       allocate(nigrid(iphi,itheta,3))

       call shf_curv_up(nterms, itheta, iphi, iup, 
     $ xmpole, ympole, zmpole, 
     $ higrid, wsigrid, nigrid)
c       write(82,*) "iup =", iup
c       do i = 1,itheta
c        do j = 1, iphi
c          write(82,*) i,j,"H",higrid(j,i),"W",wigrid(j,i) 
c        enddo
c      enddo
       allocate(ctheta(itheta), whts(itheta), theta(itheta))
       allocate(ynms(0:iterms,0:iterms,itheta/2+1))
       allocate(dwsave(4*iphi+15))
      
      !! load ctheta, whts, ynms from up_...
      do i = 1, itheta
        ctheta(i) = up_ctheta(i,iup) 
        whts(i) = up_whts(i,iup)
        theta(i) = up_theta(i,iup)
      enddo
      do j = 1, iphi*4+15
        dwsave(j) = up_dwsave(j,iup)
      enddo
      do n = 0, iterms
        do m = 0, n
          do i = 1, itheta/2+1
             ynms(n,m,i) = up_ynms(n,m,i,iup)
          enddo
        enddo
      enddo 
      
      call sphtrans_fwd_real_rsv(iterms,himpole,iphi,itheta,higrid,
     $ ctheta,whts,ynms,dwsave) 

      call eval_ups_tol_real(iterms, nterms, hmpole, himpole, res)

      iwrite = idrop * 100 + 84
      write(iwrite,*) "istep=",istep,"iup=",iup,
     $ "iterms=",iterms,"differRatioForNterms=",res

      !! trunc to order nterms hmpole
      do n = 0,nterms
        do m = 0,n
          hmpole(n, m) = himpole(n, m)    
        enddo
      enddo  

      !! if res < tol do the follow then break
      
      if (iup .eq. maxup_crv) then

        do i = 1, itheta
         do j = 1, iphi

          wsigrid(j,i) = wsigrid(j,i) / dsin(theta(i))

         enddo
        enddo

        call sphtrans_fwd_real(iterms,wsimpole,iphi,itheta,wsigrid,
     $ ctheta,whts,ynms,dwsave)   

       do k = 1, 3
        call sphtrans_fwd_real(iterms,nimpole(:,:,k),iphi,itheta,
     $ nigrid(:,:,k),ctheta,whts,ynms,dwsave)   
       enddo

         !! trunc to order nterms wmpole, nmpole
        do n = 0,nterms
          do m = 0,n
            wsmpole(n, m) = wsimpole(n, m)    
            do k = 1, 3
              nmpole(n, m, k) = nimpole(n, m, k)   
            enddo
          enddo
        enddo  

      endif
  
      deallocate(ctheta, theta, whts, ynms, dwsave)
      deallocate(himpole, wsimpole, nimpole)
      deallocate(higrid, wsigrid, nigrid)

      enddo

      !! evaluate the geometry properties at the original ntheta grid
      !! use the truncated expansion (i.e. original order)
      allocate(ctheta(ntheta), whts(ntheta), theta(ntheta))
      allocate(ynms(0:nterms,0:nterms,ntheta/2+1))
      allocate(dwsave(4*nphi+15))

      iup = 1  !! dont forget reset iup
      do i = 1, ntheta 
        ctheta(i) = up_ctheta(i,iup) 
        whts(i) = up_whts(i,iup)
        theta(i) = up_theta(i,iup)
      enddo
      do j = 1, nphi*4+15
        dwsave(j) = up_dwsave(j,iup)
      enddo
      do n = 0, nterms
        do m = 0, n
          do i = 1, ntheta/2+1
             ynms(n,m,i) = up_ynms(n,m,i,iup)
          enddo
        enddo
      enddo 

      call sphtrans_real(nterms,hmpole,nphi,ntheta,hgrid,
     $ ctheta,ynms,dwsave)    

      call sphtrans_real(nterms,wsmpole,nphi,ntheta,wsgrid,
     $ ctheta,ynms,dwsave)     

      do k = 1, 3
      call sphtrans_real(nterms,nmpole(:,:,k),nphi,ntheta,ngrid(:,:,k),
     $ ctheta,ynms,dwsave)     
      enddo

      deallocate(ctheta, whts, ynms, dwsave, theta)
      END SUBROUTINE
