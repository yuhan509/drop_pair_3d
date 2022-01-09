      subroutine nearSing_upsamp_prep(iup, iterms, itheta, iphi,
     $  xmpole, ympole, zmpole,
     $  xigrid, yigrid, zigrid, wsigrid, higrid, nigrid)
      use mod_upsamp
      implicit real *8 (a-h,o-z)

      integer :: iterms, itheta, iphi, iup
      complex *16 :: xmpole(0:nterms,0:nterms),
     $ ympole(0:nterms,0:nterms), zmpole(0:nterms,0:nterms)
      real *8 :: xigrid(iphi,itheta),yigrid(iphi,itheta),
     $ zgrid(iphi,itheta), wsigrid(iphi,itheta)
      real *8 :: higrid(iphi,itheta), nigrid(iphi,itheta,3)
 
      complex *16,dimension(:,:),allocatable :: ximpole,yimpole,zimpole
      real* 8, ALLOCATABLE:: ctheta(:), theta(:), whts(:), ynms(:,:,:) 
      complex *16, ALLOCATABLE:: dwsave(:)    
      

!       write(*,*) "nterms =", iterms, "ntheta =", itheta, "nphi =",iphi  
       allocate(ximpole(0:iterms,0:iterms))
       allocate(yimpole(0:iterms,0:iterms))
       allocate(zimpole(0:iterms,0:iterms))

       !! zero padding x,y,zmpole
       call upsamp_zero_pad_real(nterms, iterms, xmpole, ympole, zmpole,
     $ ximpole, yimpole, zimpole) 

       call shf_curv_up(nterms, itheta, iphi, iup, 
     $ xmpole, ympole, zmpole, 
     $ higrid, wsigrid, nigrid)



       allocate(ctheta(itheta), theta(itheta), whts(itheta))
       allocate(ynms(0:iterms,0:iterms,itheta/2+1))
       allocate(dwsave(4*iphi+15))
      
      !! load ctheta, whts, ynms from up_...
      do i = 1, itheta
        ctheta(i) = up_ctheta(i,iup) 
        whts(i) = up_whts(i,iup)
        theta(i) = up_theta(i,iup) !! needed in test
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
      do i = 1, itheta
        do j = 1, iphi
            wsigrid(j,i) = wsigrid(j,i) / dsin(theta(i))
        enddo
      enddo
      !! transform x,y,zimpole to xigrid, yigrid, zigrid
      call sphtrans_real(iterms,ximpole,iphi,itheta,  
     $  xigrid,ctheta,ynms,dwsave)
      call sphtrans_real(iterms,yimpole,iphi,itheta,  
     $  yigrid,ctheta,ynms,dwsave)
      call sphtrans_real(iterms,zimpole,iphi,itheta,  
     $  zigrid,ctheta,ynms,dwsave)   

       deallocate(ctheta, theta, whts, ynms, dwsave)
       deallocate(ximpole,yimpole,zimpole)
      END SUBROUTINE
c
c
c
c
      subroutine nearSing_upsamp_En(iup, iterms, itheta, iphi,
     $  en1mpole, en2mpole, En1igrid, En2igrid)
      use mod_upsamp
      implicit real *8 (a-h,o-z)

      integer :: iterms, itheta, iphi, iup
      complex *16 :: en1mpole(0:nterms,0:nterms),
     $ en2mpole(0:nterms,0:nterms)
      real *8 :: en1igrid(iphi,itheta),en2igrid(iphi,itheta)

      complex *16,dimension(:,:),allocatable :: en1impole,en2impole
      real* 8, ALLOCATABLE:: ctheta(:), theta(:), whts(:), ynms(:,:,:) 
      complex *16, ALLOCATABLE:: dwsave(:)    
      
!       write(*,*) "nterms =", iterms, "ntheta =", itheta, "nphi =",iphi  
       allocate(en1impole(0:iterms,0:iterms))
       allocate(en2impole(0:iterms,0:iterms))

       !! zero padding 
       do n = 0, nterms
        do m = 0, n
           en1impole(n,m) = en1mpole(n,m)
           en2impole(n,m) = en2mpole(n,m)
        enddo
       enddo
       do n = nterms+1, iterms
        do m = 0, n
           en1impole(n,m) = (0d0,0d0)
           en2impole(n,m) = (0d0,0d0)
        enddo
       enddo       

       allocate(ctheta(itheta), theta(itheta), whts(itheta))
       allocate(ynms(0:iterms,0:iterms,itheta/2+1))
       allocate(dwsave(4*iphi+15))
      
      !! load ctheta, whts, ynms from up_...
      do i = 1, itheta
        ctheta(i) = up_ctheta(i,iup) 
        whts(i) = up_whts(i,iup)
        theta(i) = up_theta(i,iup) !! needed in test
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

      !! transform x,y,zimpole to xigrid, yigrid, zigrid
      call sphtrans_real(iterms,en1impole,iphi,itheta,  
     $  en1igrid,ctheta,ynms,dwsave)
      call sphtrans_real(iterms,en2impole,iphi,itheta,  
     $  en2igrid,ctheta,ynms,dwsave)

       deallocate(ctheta, theta, whts, ynms, dwsave)
       deallocate(en1impole, en2impole)
      END SUBROUTINE
c
c
c
      subroutine nearSing_upsamp_prep_mpole(iup, iterms, itheta, iphi,
     $  xmpole, ympole, zmpole, hmpole, wsmpole, nmpole, 
     $  xigrid, yigrid, zigrid, wsigrid, higrid, nigrid)
      use mod_upsamp
      implicit real *8 (a-h,o-z)

      integer :: iterms, itheta, iphi, iup
      complex *16 :: xmpole(0:nterms,0:nterms),
     $ ympole(0:nterms,0:nterms), zmpole(0:nterms,0:nterms)
      complex *16 :: hmpole(0:nterms,0:nterms),
     $ wsmpole(0:nterms,0:nterms), nmpole(0:nterms,0:nterms,3)
      real *8 :: xigrid(iphi,itheta),yigrid(iphi,itheta),
     $ zgrid(iphi,itheta), wsigrid(iphi,itheta)
      real *8 :: higrid(iphi,itheta), nigrid(iphi,itheta,3)
 
      complex *16,dimension(:,:),allocatable :: ximpole,yimpole,zimpole
      complex *16,dimension(:,:),allocatable :: himpole,wsimpole
      complex *16,dimension(:,:,:),allocatable :: nimpole
      real* 8, ALLOCATABLE:: ctheta(:), theta(:), whts(:), ynms(:,:,:) 
      complex *16, ALLOCATABLE:: dwsave(:)    
      

!       write(*,*) "nterms =", iterms, "ntheta =", itheta, "nphi =",iphi  
       allocate(ximpole(0:iterms,0:iterms))
       allocate(yimpole(0:iterms,0:iterms))
       allocate(zimpole(0:iterms,0:iterms))
       allocate(himpole(0:iterms,0:iterms))
       allocate(wsimpole(0:iterms,0:iterms))
       allocate(nimpole(0:iterms,0:iterms,3))

       !! zero padding x,y,zmpole
       call upsamp_zero_pad_real(nterms, iterms, xmpole, ympole, zmpole,
     $ ximpole, yimpole, zimpole) 

       !! zero padding h,w,mpole
       do n = 0, iterms
         do m = 0, n
           himpole(n,m) = (0d0, 0d0)
           wsimpole(n,m) = (0d0, 0d0)
           do k = 1, 3
             nimpole(n,m,k) = (0d0, 0d0)
           enddo
         enddo
       enddo
       do n = 0, nterms
        do m = 0, n
          himpole(n,m) = hmpole(n,m)
          wsimpole(n,m) = wsmpole(n,m)
          do k = 1, 3
            nimpole(n,m,k) = nmpole(n,m,k)
          enddo          
        enddo
       enddo


       allocate(ctheta(itheta), theta(itheta), whts(itheta))
       allocate(ynms(0:iterms,0:iterms,itheta/2+1))
       allocate(dwsave(4*iphi+15))
      
      !! load ctheta, whts, ynms from up_...
      do i = 1, itheta
        ctheta(i) = up_ctheta(i,iup) 
        whts(i) = up_whts(i,iup)
        theta(i) = up_theta(i,iup) !! needed in test
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

      !! transform x,y,zimpole to xigrid, yigrid, zigrid
      call sphtrans_real(iterms,ximpole,iphi,itheta,  
     $  xigrid,ctheta,ynms,dwsave)
      call sphtrans_real(iterms,yimpole,iphi,itheta,  
     $  yigrid,ctheta,ynms,dwsave)
      call sphtrans_real(iterms,zimpole,iphi,itheta,  
     $  zigrid,ctheta,ynms,dwsave)   


      call sphtrans_real(iterms,himpole,iphi,itheta,  
     $  higrid,ctheta,ynms,dwsave)
      call sphtrans_real(iterms,wsimpole,iphi,itheta,  
     $  wsigrid,ctheta,ynms,dwsave)
      
      do k = 1, 3
        call sphtrans_real(iterms,nimpole(:,:,k),iphi,itheta,  
     $  nigrid(:,:,k),ctheta,ynms,dwsave)
      enddo

    !  do i = 1, itheta
    !    do j = 1, iphi
    !       wsigrid(j,i) = wsigrid(j,i) / dsin(theta(i))
    !    enddo
    !  enddo

 !     write(9102,*) "j,i,hgrid,wsgrid,ngrid"  
 !     do i = 1, itheta
 !      do j = 1, iphi
 !       write(9102,*) j,i,higrid(j,i),wsigrid(j,i),nigrid(j,i,1)  
 !    $ ,nigrid(j,i,2),nigrid(j,i,3)/ctheta(i)
 !       write(9102,*) nigrid(j,i,1)**2+nigrid(j,i,2)**2+nigrid(j,i,3)**2
 !      enddo
 !     enddo

       deallocate(ctheta, theta, whts, ynms, dwsave)
       deallocate(ximpole,yimpole,zimpole)
      END SUBROUTINE