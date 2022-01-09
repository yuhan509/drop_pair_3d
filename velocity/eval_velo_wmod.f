      SUBROUTINE eval_velo(istep, Ca,
     $ xgrid, ygrid, zgrid, xmpole, ympole, zmpole,
     $ hgrid, wgrid, ngrid,
     $ ugrids, enmpole)
     
      use mod_orgHelper
      use mod_rotsphf
      
      implicit real *8 (a-h,o-z)
      
      integer :: istat = 2, istep
      real *8 :: Ca
      real *8 :: xgrid(nphi,ntheta),ygrid(nphi,ntheta),
     $ zgrid(nphi,ntheta)

      complex *16 :: xmpole(0:nterms,0:nterms),
     $ ympole(0:nterms,0:nterms),zmpole(0:nterms,0:nterms)
      real *8 :: hgrid(nphi,ntheta), wgrid(nphi,ntheta), 
     $ ngrid(nphi,ntheta,3)
      !! output 
      real *8 :: ugrids(nphi,ntheta,3)
      complex *16 :: enmpole(0:nterms,-nterms:nterms)

      real* 8, ALLOCATABLE:: swhts(:)
      complex* 16 :: ctmp,cquadr,ctmp1,ctmp2,ctmp3,
     $ cquadr1,cquadr2,cquadr3
      real *8, ALLOCATABLE:: Einfngrid(:,:)

      complex *16, ALLOCATABLE:: mpole(:,:), fgrid(:,:), ggrid(:,:)
      complex *16, pointer, dimension(:,:,:,:):: grids
      real *8, pointer, dimension(:,:,:,:):: rxgrids,rygrids,rzgrids
      complex *16, pointer, dimension(:,:,:,:):: clapl_quadrs
      complex *16, pointer, dimension(:,:,:,:):: coef_mpole
      complex *16, ALLOCATABLE:: cap_lpl(:,:),cap_Einf(:)
      complex *16, ALLOCATABLE:: testa(:,:),testb(:,:)
      INTEGER :: NRHS=1,N2
      INTEGER :: INFO
      INTEGER,allocatable :: ipiv(:),ipivtest(:)
      
      complex *16 :: eiphi,conj_eiphi
      real*8, ALLOCATABLE:: En(:,:)
      real *8, ALLOCATABLE :: fx(:,:),fy(:,:),fz(:,:),ws(:,:)
      real *8, ALLOCATABLE :: fxgrids(:,:,:,:), fygrids(:,:,:,:),
     $ fzgrids(:,:,:,:),wsgrids(:,:,:,:)
      complex *16, ALLOCATABLE :: fxMpole(:,:),
     $ fyMpole(:,:),fzMpole(:,:), wsmpole(:,:)
      real *8, ALLOCATABLE :: stk_quadrs(:,:,:)
      
      N2=(nterms+1)*(nterms+1)

      allocate(swhts(ntheta))
      allocate(mpole(0:nterms,-nterms:nterms))
      allocate(fgrid(nphi,ntheta))
      allocate(ggrid(nphi,ntheta))
ccc   rotation related:
      allocate(grids(nphi,ntheta,nphi,ntheta))
      allocate(rxgrids(nphi,ntheta,nphi,ntheta))      
      allocate(rygrids(nphi,ntheta,nphi,ntheta))
      allocate(rzgrids(nphi,ntheta,nphi,ntheta))
      allocate(Einfngrid(nphi,ntheta))
      allocate(clapl_quadrs(nphi,ntheta,0:nterms,-nterms:nterms))
      allocate(coef_mpole(0:nterms,-nterms:nterms,
     $ 0:nterms,-nterms:nterms))
!      allocate(cap_lpl(N2-1,N2-1))
!      allocate(cap_Einf(N2-1))
!      allocate(ipiv(N2-1))
      allocate(cap_lpl(N2,N2))   !! n start with 0
      allocate(cap_Einf(N2))     !! n start with 0 
      allocate(ipiv(N2))         !! n start with 0
      
      allocate(En(nphi,ntheta))
      allocate(fx(nphi,ntheta))
      allocate(fy(nphi,ntheta))
      allocate(fz(nphi,ntheta))
      allocate(ws(nphi,ntheta))
      allocate(fxgrids(nphi,ntheta,nphi,ntheta))
      allocate(fygrids(nphi,ntheta,nphi,ntheta))
      allocate(fzgrids(nphi,ntheta,nphi,ntheta))
      allocate(wsgrids(nphi,ntheta,nphi,ntheta))      
      allocate(fxMpole(0:nterms,0:nterms))
      allocate(fyMpole(0:nterms,0:nterms))
      allocate(fzMpole(0:nterms,0:nterms))     
      allocate(wsmpole(0:nterms,0:nterms))   
      allocate(stk_quadrs(nphi,ntheta,3))

!      write(*,*) nphi,ntheta,'ynms',ynms(3,2,1), ctheta(3), zwsave(10)
!      write(*,*) sphfgrid(1,2,3,3),rsphfgrids(1,2,2,2,2,2)           

      PI = 4.D0*DATAN(1.D0)
cccccccccccccccccccccccccccccccccccccccccccccc
      do i=1,ntheta
        do j=1,nphi
          !! E_inf * ngrid
          ws(j,i) = wgrid(j,i)/dsin(theta(i))
          Einfngrid(j,i) = ngrid(j,i,3)
          write(88,*) j, i, ngrid(j,i,3), ctheta(i)
        enddo
      enddo

      call sphtrans_fwd_real_rsv(nterms,wsmpole,nphi,ntheta,ws,
     $  ctheta,whts,ynms,dwsave)
cccccccccccccccccccccccccccccccccccc
cccccc rotate x,y,z mpole and local area element ccccccccccc
      nbeta = ntheta
      nrot = nphi
      call rotgrid_dsr_real(nterms,xmpole,nphi,ntheta,
     $ nbeta,theta,nrot,rxgrids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,ympole,nphi,ntheta,
     $ nbeta,theta,nrot,rygrids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,zmpole,nphi,ntheta,
     $ nbeta,theta,nrot,rzgrids,
     $ rotmat,ctheta,ynms,dwsave)     
      call rotgrid_dsr_real(nterms,wsmpole,nphi,ntheta,
     $ nbeta,theta,nrot,wsgrids,
     $ rotmat,ctheta,ynms,dwsave)     

cccccccccccccccccccccccccccccccccccccccccccccc
cc    special gauss-legendre quadrature weights
cccccccccccccccccccccccccccccccccccccccccccccc      
      do i = 1, ntheta
        swhts(i) = 0d0
        do n = 0, nterms
          call legepol(ctheta(i),n,pol,der)
c          write(*,*) 'x=',ctheta(j),'n=',i,'Pn(x)=',pol
          swhts(i) = swhts(i) + pol !/dsqrt(4d0*i + 2d0)
        enddo          
        swhts(i) = swhts(i) * whts(i) * 2d0 * dsin(theta(i)/2d0)
      enddo    

cc    special gauss-legendre quadrature(lapalce) for all rotated grids
      
      do n = 0,nterms
        do m = -n,n
          !!write(*,*) "using rotating sph", n, m
          do kk = 1,nrot
            do k = 1,nbeta
              cquadr = (0d0,0d0)
              do i = 1,ntheta
                ctmp = (0d0,0d0)
                do j = 1, nphi
                  dist = (xgrid(kk,k)-rxgrids(j,i,kk,k))**2 
     $             + (ygrid(kk,k)-rygrids(j,i,kk,k))**2
     $             + (zgrid(kk,k)-rzgrids(j,i,kk,k))**2
                  dist = dsqrt(dist)
              tmp_inner = ngrid(kk,k,1)*(xgrid(kk,k)-rxgrids(j,i,kk,k))
     $ + ngrid(kk,k,2)*(ygrid(kk,k)-rygrids(j,i,kk,k)) 
     $ + ngrid(kk,k,3)*(zgrid(kk,k)-rzgrids(j,i,kk,k)) 
              ctmp = ctmp + tmp_inner/(dist**3)*rsphfgrids(j,i,kk,k,n,m)
     $ * wsgrids(j,i,kk,k)
!                write(*,*) "tmp_inner/(dist**3)",tmp_inner/(dist**3)
!                write(*,*) "1/dist", 1/dist
                enddo
                cquadr = ctmp * swhts(i) + cquadr
              enddo    
              clapl_quadrs(kk,k,n,m) = cquadr/(1d0*nphi)*2*pi/(-4*pi)
!              write(13,*) 'n=',n,'m=',m,'irot=',kk,'jbeta=',k,
!     $  'special g-l quadr laplace',clapl_quadrs(kk,k,n,m)
            enddo
          enddo
        enddo
      enddo     

ccccccccc Inner Product ccccccccccccc

      do n = 0,nterms
        do m = -n,n
         do n1 = 0,nterms
         do m1 = -n1, n1
           cquadr = (0d0,0d0)
           do k = 1, ntheta
             ctmp = (0d0,0d0)
             do l = 1, nphi
               ctmp = ctmp + dconjg( sphfgrid(l,k,n1,m1) )
     $ * clapl_quadrs(l,k,n,m)
             enddo
             cquadr = cquadr + ctmp * whts(k)
           enddo
           coef_mpole(n1,m1,n,m) = cquadr/dble(nphi)*2d0*pi
!              write(14,*) 'n=',n,'m=',m,'n1=',n1,'m1=',m1,
!     $ coef_mpole(n1,m1,n,m)
         enddo
       enddo     
      enddo
      enddo  

       do n = 0,nterms
         do m = -n, n
           cquadr = (0d0,0d0)
           do k = 1, ntheta
             ctmp = (0d0,0d0)
             do l = 1, nphi
            ctmp=ctmp+dconjg( sphfgrid(l,k,n,m) )*Einfngrid(l,k)       
             enddo
             cquadr = cquadr + ctmp * whts(k)
           enddo
           mpole(n, m) = cquadr/dble(nphi)*2d0*pi
         enddo
       enddo
       
       
!!      do n = 0,nterms
!!        do m = -n,n
!!         do n1 = 0,nterms
!!         do m1 = -n1, n1
!!           cquadr = (0d0,0d0)
!!           do k = 1, ntheta
!!             ctmp = (0d0,0d0)
!!             do l = 1, nphi
!!               ctmp = ctmp + dconjg( sphfgrid(l,k,n1,m1) )
!!     $ * sphfgrid(l,k,n,m)
!!             enddo
!!             cquadr = cquadr + ctmp * whts(k)
!!           enddo
!!           write(*,*) n,m,n1,m1,cquadr/dble(nphi)*2d0*pi
!!!              write(14,*) 'n=',n,'m=',m,'n1=',n1,'m1=',m1,
!!!     $ coef_mpole(n1,m1,n,m)
!!         enddo
!!       enddo     
!!      enddo
!!      enddo  
cccccccccccccccccccccccccccccccccccc
      j = 1
      do n = 0,nterms !! start 0
        do m = -n, n
          i = 1
          do n1 = 0,nterms !! start 0
            do m1 = -n1, n1
              cap_lpl(i,j)=coef_mpole(n1,m1,n,m)
!              write(15,*) 'n=',n,'m=',m,'n1=',n1,'m1=',m1,'i=',i,'j=',j,
!     $  cap_lpl(i,j)
              i = i + 1    
            enddo
          enddo
          j = j + 1
        enddo
      enddo    
      
      
    
      do i = 1, j-1
          cap_lpl(i,i) = cap_lpl(i,i)+dcmplx(0.5d0*4*pi,0d0)
!         write(15, *) 'i=',i,cap_lpl(i,i)
      enddo

      write(14,*) 'Y proj of E_inf*ngrid, ALL B terms in A * X = B'
      write(21,*)  "only nonzero B terms, in A * X = B"
      i = 1
      do n = 0, nterms !! start 0
          do m = -n, n              
              cap_Einf(i) = mpole(n, m)
             ! write(14,*) 'n',n,'m',m,'i=',i,cap_Einf(i)
              if (abs(cap_Einf(i)*dconjg(cap_Einf(i))).lt.1d-25) then
                cap_Einf(i) = 0d0
              else
                 write(21,*)  'n',n,'m',m,'i=',i,cap_Einf(i)
              endif
              i = i + 1
          enddo
      enddo
      
cccccccccccc
       write(19,*)  " all A terms, in A * X = B"   
       write(20,*)  " only nonzero A terms, in A * X = B"  
       do k = 1, N2
         do l = 1, N2
           !write(19,*) 'i=',k,'j=',l, cap_lpl(k,l)
           if (abs(cap_lpl(k,l)*dconjg(cap_lpl(k,l))).lt.1d-25) then 
             cap_lpl(k,l) = 0d0
           else
            write(20,*) 'i=',k,'j=',l, cap_lpl(k,l)
           endif
          enddo
       enddo
     
cccccccccccc
            
      CALL zgesv(N2,NRHS,cap_lpl,N2,IPIV,cap_Einf,N2,INFO)      
!      CALL zgesv(N2-1,NRHS,cap_lpl,N2-1,
!     $ IPIV,cap_Einf,N2-1,INFO)      
      write(16,*) 'linear solver zgesv solution:','INFO=',INFO
      
      i = 1
      do n = 0, nterms !! start 0
          do m = -n, n              
              enmpole(n, m) = cap_Einf(i)   
              write(16,*) 'i=',i,'n',n,'m',m,cap_Einf(i)
              if (abs(cap_Einf(i)*dconjg(cap_Einf(i))).lt.1d-25) then
                enmpole(n, m) = 0d0
              else
              
              endif
              i = i + 1
          enddo
      enddo
      
ccccc spherical harmonic transform ccccccc
      !!!! special case for the first step sphere shape
      if (istep .eq. 1) then 
        enmpole(0,0) = (0d0, 0d0)
      endif

      do j = 1,ntheta
          do i = 1, nphi
            ctmp = (0d0,0d0)
            do n = 0, nterms   !! start 0, only take odd terms
              do m = -n, n              
                  ctmp = ctmp + enmpole(n, m)*sphfgrid(i,j,n,m)             
              enddo
            enddo 
            En(i,j) = dble(ctmp)
            write(17,*) 'j=',j,'i=',i,En(i,j),En(i,j)/dcos(theta(j))
          enddo
      enddo

cccccccccc total force mpole cccccccccc
c          ( En^2/2*Ca + 2*H ) * ngrid
ccccccccccccccccccccccccccccccccccccccc    

      Caflow = 0.0d0
      do i = 1, ntheta
        do j = 1, nphi            
          tmp = (En(j,i)**2/2d0*Ca + 2d0*hgrid(j,i))
          fx(j,i) = tmp*ngrid(j,i,1)
          fy(j,i) = tmp*ngrid(j,i,2)
          fz(j,i) = tmp*ngrid(j,i,3)
          write(32,*) j,i,En(j,i)/dcos(theta(i)),dcos(theta(i)),ws(j,i)
        enddo
      enddo
      
      call sphtrans_fwd_real(nterms,fxmpole,nphi,ntheta,fx,
     $     ctheta,whts,ynms,dwsave)
      call sphtrans_fwd_real(nterms,fympole,nphi,ntheta,fy,
     $     ctheta,whts,ynms,dwsave)    
      call sphtrans_fwd_real(nterms,fzmpole,nphi,ntheta,fz,
     $     ctheta,whts,ynms,dwsave)           

cccccccccccccccccccccccccccccccccccccccc

cccccc rotate force components mpole ccccccccccc
      call rotgrid_dsr_real(nterms,fxmpole,nphi,ntheta,
     $ nbeta,theta,nrot,fxgrids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,fympole,nphi,ntheta,
     $ nbeta,theta,nrot,fygrids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,fzmpole,nphi,ntheta,
     $ nbeta,theta,nrot,fzgrids,
     $ rotmat,ctheta,ynms,dwsave)            

cccccccccccccccccccccccccccccccccccccccccccccc

cc    special gauss-legendre quadrature(lapalce) for all rotated grids
          write(34,*) 'irot ','jbeta ', 'special g-l quadr stk'
          do kk = 1,nrot
            do k = 1,nbeta
              quadr1 = 0d0
              quadr2 = 0d0
              quadr3 = 0d0        
              do i = 1,ntheta
                tmp1 = 0d0
                tmp2 = 0d0
                tmp3 = 0d0
                do j = 1, nphi
                  dist = (xgrid(kk,k)-rxgrids(j,i,kk,k))**2 
     $             + (ygrid(kk,k)-rygrids(j,i,kk,k))**2
     $             + (zgrid(kk,k)-rzgrids(j,i,kk,k))**2
                  dist = dsqrt(dist)
              tmp = (xgrid(kk,k)-rxgrids(j,i,kk,k)) * fxgrids(j,i,kk,k)
     $ + (ygrid(kk,k)-rygrids(j,i,kk,k))* fygrids(j,i,kk,k)
     $ + (zgrid(kk,k)-rzgrids(j,i,kk,k))* fzgrids(j,i,kk,k)
              tmp1 = tmp1 + (1d0/dist * fxgrids(j,i,kk,k)
     $ + (xgrid(kk,k)-rxgrids(j,i,kk,k))/dist**3 * tmp) 
     $ *wsgrids(j,i,kk,k)
              tmp2 = tmp2 + (1d0/dist * fygrids(j,i,kk,k)
     $ + (ygrid(kk,k)-rygrids(j,i,kk,k))/dist**3 * tmp)
     $ *wsgrids(j,i,kk,k) 
              tmp3 = tmp3 + (1d0/dist * fzgrids(j,i,kk,k)
     $ + (zgrid(kk,k)-rzgrids(j,i,kk,k))/dist**3 * tmp) 
     $ *wsgrids(j,i,kk,k)  
                enddo
                quadr1 = tmp1 * swhts(i) + quadr1
                quadr2 = tmp2 * swhts(i) + quadr2
                quadr3 = tmp3 * swhts(i) + quadr3
              enddo    
              !>! weight for quadrature over phi is 2*pi/nphi combine
              !>! with factor 1/(8*pi) for evaluating the velocity for
              !>! a viscosity ratio one
              ugrids(kk,k,1) = quadr1/(4*nphi) - Caflow*xgrid(kk,k)/2d0
              ugrids(kk,k,2) = quadr2/(4*nphi) - Caflow*ygrid(kk,k)/2d0
              ugrids(kk,k,3) = quadr3/(4*nphi) + Caflow*zgrid(kk,k)
              write(34,*) kk, k, ugrids(kk,k,1), ugrids(kk,k,2),
     $  ugrids(kk,k,3)
            enddo
          enddo 
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




