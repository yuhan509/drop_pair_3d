      subroutine upsamp_crv_pre_cal()
     
      use mod_upsamp
     
      implicit real *8 (a-h,o-z)
      integer :: iterms, iphi, itheta

      complex *16, pointer, dimension(:,:,:,:) :: yfgrid,dydth, 
     $ dydph,d2ydth2,d2ydthdph,d2ydph2
      real* 8, ALLOCATABLE:: ctheta(:), whts(:), ynms(:,:,:), theta(:)
      complex *16, ALLOCATABLE::  zwsave(:), dwsave(:), mpole(:,:)

       do iup = 1, maxup_crv
        iterms = iup * nterms
        call grid_dim_init(iterms,itheta,iphi)

        allocate(ctheta(itheta), whts(itheta), theta(itheta))
        allocate(ynms(0:iterms,0:iterms,itheta/2+1))
        allocate(zwsave(4*iphi+15),dwsave(4*iphi+15))
        allocate(mpole(0:iterms, -iterms:iterms))
        
        allocate(yfgrid(iphi,itheta,0:nterms,-nterms:nterms))
        allocate(dydth(iphi,itheta,0:nterms,-nterms:nterms))
        allocate(dydph(iphi,itheta,0:nterms,-nterms:nterms))
        allocate(d2ydth2(iphi,itheta,0:nterms,-nterms:nterms))
        allocate(d2ydthdph(iphi,itheta,0:nterms,-nterms:nterms))
        allocate(d2ydph2(iphi,itheta,0:nterms,-nterms:nterms))
        
        do n = 0, iterms
         do m = -n, n      
          mpole(n,m) = (0d0, 0d0)
         enddo
        enddo
         write(*,*) "ccc1",iterms
         call grid_lege_init(iterms,itheta,iphi,
     $    ctheta, theta, whts, ynms, dwsave, zwsave)

         write(*,*) "ccc2",iterms
c         call grid_sphf(iterms, itheta, iphi, 
c     $   ctheta, ynms, zwsave, yfgrid) 

        do n = 0, nterms
         do m = -n, n
          mpole(n,m) = (1d0, 0d0)
          call sphtrans_cmpl(iterms,mpole,iphi,itheta,yfgrid(:,:,n,m),
     $     ctheta,ynms,zwsave)          
          mpole(n,m) = (0d0, 0d0)
         enddo
        enddo
          write(*,*) "ccc3",iterms
         call shf_diffs(yfgrid, nterms, iphi, itheta, 
     $ dydth, dydph, d2ydth2, d2ydthdph, d2ydph2)
          write(*,*) "ccc4",iterms

          do n = 0, nterms
            do m = 0, n
              do i = 1, itheta
                do j = 1, iphi
                   ! write(*,*) iup,j,i,n,m
                  up_sphfgrid(j,i,n,m,iup) = yfgrid(j,i,n,m)
                  up_dydph(j,i,n,m,iup) = dydph(j,i,n,m)
                  up_dydth(j,i,n,m,iup) = dydth(j,i,n,m)
                  up_d2ydph2(j,i,n,m,iup) = d2ydph2(j,i,n,m)
                  up_d2ydth2(j,i,n,m,iup) = d2ydth2(j,i,n,m)
                  up_d2ydthdph(j,i,n,m,iup) = d2ydthdph(j,i,n,m)                 
                enddo
              enddo   
            enddo
          enddo 

          write(*,*) "ccc5",iterms
          do i = 1, itheta
            up_ctheta(i,iup) = ctheta(i)
            up_theta(i,iup) = theta(i)
            up_whts(i,iup) = whts(i)
          enddo
          write(*,*) "ccc6",iterms
          do j = 1, iphi*4+15
            up_zwsave(j,iup) = zwsave(j)
            up_dwsave(j,iup) = dwsave(j)
          enddo
          write(*,*) "ccc7",iterms
          do n = 0, iterms
            do m = 0, iterms
              do i = 1, itheta/2+1
                up_ynms(n,m,i,iup) = ynms(n,m,i)
              enddo
            enddo
          enddo 
          write(*,*) "ccc8",iterms
          deallocate(ctheta, whts, theta, ynms, zwsave, dwsave, mpole)
          deallocate(yfgrid,dydth,dydph,d2ydth2,d2ydthdph,d2ydph2)
       enddo       
       
      END subroutine
c
c
c
c
      subroutine upsamp_upx_pre_cal()
     
      use mod_upsamp
     
      implicit real *8 (a-h,o-z)
      integer :: iterms, iphi, itheta

      real* 8, ALLOCATABLE:: ctheta(:), whts(:), ynms(:,:,:), theta(:)
      complex *16, ALLOCATABLE::  zwsave(:), dwsave(:), mpole(:,:)

        iup = maxup 
        iterms = iup * nterms
        call grid_dim_init(iterms,itheta,iphi)

        allocate(ctheta(itheta), whts(itheta), theta(itheta))
        allocate(ynms(0:iterms,0:iterms,itheta/2+1))
        allocate(zwsave(4*iphi+15),dwsave(4*iphi+15))
        allocate(mpole(0:iterms, -iterms:iterms))
        
        do n = 0, iterms
         do m = -n, n      
          mpole(n,m) = (0d0, 0d0)
         enddo
        enddo
         write(*,*) "ccc1",iterms
         call grid_lege_init(iterms,itheta,iphi,
     $    ctheta, theta, whts, ynms, dwsave, zwsave)

         write(*,*) "ccc2",iterms

        do n = 0, nterms
         do m = 0, n
          mpole(n,m) = (1d0, 0d0)
          call sphtrans_cmpl(iterms,mpole,iphi,itheta,
     $     upx_sphfgrid(:,:,n,m),ctheta,ynms,zwsave)   
          write(*,*) "build upx_sphfgrid(:,:,n,m)",n,m      
          mpole(n,m) = (0d0, 0d0)
         enddo
        enddo
          write(*,*) "ccc3",iterms

          write(*,*) "ccc5",iterms
          do i = 1, itheta
            up_ctheta(i,iup) = ctheta(i)
            up_theta(i,iup) = theta(i)
            up_whts(i,iup) = whts(i)
          enddo
          write(*,*) "ccc6",iterms
          do j = 1, iphi*4+15
            up_zwsave(j,iup) = zwsave(j)
            up_dwsave(j,iup) = dwsave(j)
          enddo
          write(*,*) "ccc7",iterms
          do n = 0, iterms
            do m = 0, iterms
              do i = 1, itheta/2+1
                up_ynms(n,m,i,iup) = ynms(n,m,i)
              enddo
            enddo
          enddo 
          write(*,*) "ccc8",iterms
          deallocate(ctheta, whts, theta, ynms, zwsave, dwsave, mpole)
    
      END subroutine
c     
c     
c
c      
      subroutine upsamp_repara_pre_cal()
     
      use mod_upsamp
     
      implicit real *8 (a-h,o-z)
      integer :: iterms, iphi, itheta

      complex *16, pointer, dimension(:,:,:,:) :: yfgrid,dydth, 
     $ dydph,d2ydth2,d2ydthdph,d2ydph2
      real* 8, ALLOCATABLE:: ctheta(:), whts(:), ynms(:,:,:), theta(:)
      complex *16, ALLOCATABLE::  zwsave(:), dwsave(:), mpole(:,:)

        iterms = nup_repara * nterms
        call grid_dim_init(iterms,itheta,iphi)

        allocate(ctheta(itheta), whts(itheta), theta(itheta))
        allocate(ynms(0:iterms,0:iterms,itheta/2+1))
        allocate(zwsave(4*iphi+15),dwsave(4*iphi+15))
        allocate(mpole(0:iterms, -iterms:iterms))
        
        allocate(yfgrid(iphi,itheta,0:iterms,-iterms:iterms))
        allocate(dydth(iphi,itheta,0:iterms,-iterms:iterms))
        allocate(dydph(iphi,itheta,0:iterms,-iterms:iterms))
        allocate(d2ydth2(iphi,itheta,0:iterms,-iterms:iterms))
        allocate(d2ydph2(iphi,itheta,0:iterms,-iterms:iterms))
        allocate(d2ydthdph(iphi,itheta,0:iterms,-iterms:iterms))
      
        do n = 0, iterms
         do m = -n, n      
          mpole(n,m) = (0d0, 0d0)
         enddo
        enddo
         write(*,*) "ccc1r",iterms
         call grid_lege_init(iterms,itheta,iphi,
     $    ctheta, theta, whts, ynms, dwsave, zwsave)

         write(*,*) "ccc2r",iterms
c         call grid_sphf(iterms, itheta, iphi, 
c     $   ctheta, ynms, zwsave, yfgrid) 

        do n = 0, iterms
         do m = -n, n
          mpole(n,m) = (1d0, 0d0)
          call sphtrans_cmpl(iterms,mpole,iphi,itheta,yfgrid(:,:,n,m),
     $     ctheta,ynms,zwsave)          
          mpole(n,m) = (0d0, 0d0)
         enddo
        enddo
          write(*,*) "ccc3r",iterms
         call shf_diffs(yfgrid, iterms, iphi, itheta, 
     $ dydth, dydph, d2ydth2, d2ydthdph, d2ydph2)
          write(*,*) "ccc4r",iterms

          do n = 0, iterms
            do m = 0, n
              do i = 1, itheta
                do j = 1, iphi
                   ! write(*,*) iup,j,i,n,m
                  upr_sphfgrid(j,i,n,m) = yfgrid(j,i,n,m)
                  upr_dydph(j,i,n,m) = dydph(j,i,n,m)
                  upr_dydth(j,i,n,m) = dydth(j,i,n,m)           
                enddo
              enddo   
            enddo
          enddo 

      END subroutine
