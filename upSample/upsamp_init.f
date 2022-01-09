      SUBROUTINE grid_dim_init(nterms, ntheta, nphi)
      integer :: nterms, next235, nphi, ntheta
      
      ntheta = nterms + 1
      nphi = 2 * nterms + 2
      call fftnext235(nphi, next235)
      nphi = next235      
      
      END SUBROUTINE      
c
c
c
      SUBROUTINE grid_lege_init(nterms, ntheta, nphi, ctheta, theta,
     $ whts, ynms, dwsave, zwsave)
      implicit real *8 (a-h,o-z)
        integer :: nterms, nphi, ntheta
      
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8 :: ctheta(ntheta), whts(ntheta), theta(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 :: dwsave(4*nphi+15), zwsave(4*nphi+15)
c
        allocate(xs(ntheta))
        allocate(ws(ntheta))
        allocate(u(ntheta,ntheta))
        allocate(v(ntheta,ntheta))
        allocate(rat1(nterms+1,nterms+1))
        allocate(rat2(nterms+1,nterms+1))
c
c      ... construct Gauss_Legendre nodes and weights on the interval [-1,1]
c
        itype=1
        call legeexps(itype,ntheta,xs,u,v,ws)
        do i=1,ntheta
c       
        ctheta(i)=-xs(i)
        whts(i)=ws(i)
        theta(i) = dacos(ctheta(i))
        enddo
c
c
        call ylgndrini(nterms,rat1,rat2)
c
        done=1
        pi=4*atan(done)
c
        call zffti(nphi,zwsave)
        call dffti(nphi,dwsave)
c
        do k=1,ntheta/2+1
        call ylgndrf(nterms,ctheta(k),ynms(0,0,k),rat1,rat2)
        enddo
c
        deallocate(xs, ws, u, v, rat1, rat2)
        return      
      END SUBROUTINE     
c
c
c
      SUBROUTINE grid_lege_init_up(nterms, ntheta, nphi, 
     $ ctheta, whts, ynms, dwsave, zwsave)
      implicit real *8 (a-h,o-z)
        integer :: nterms, nphi, ntheta
      
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8 :: ctheta(ntheta), whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 dwsave(4*nphi+15), zwsave(4*nphi+15)
c
        allocate(xs(ntheta))
        allocate(ws(ntheta))
        allocate(u(ntheta,ntheta))
        allocate(v(ntheta,ntheta))
        allocate(rat1(nterms+1,nterms+1))
        allocate(rat2(nterms+1,nterms+1))
c
c      ... construct Gauss_Legendre nodes and weights on the interval [-1,1]
c
        itype=1
        call legeexps(itype,ntheta,xs,u,v,ws)
        do i=1,ntheta
c       
        ctheta(i)=-xs(i)
        whts(i)=ws(i)
        enddo
c
c
        call ylgndrini(nterms,rat1,rat2)
c
        done=1
        pi=4*atan(done)
c
        call zffti(nphi,zwsave)
        call dffti(nphi,dwsave)
c
        do k=1,ntheta/2+1
        call ylgndrf(nterms,ctheta(k),ynms(0,0,k),rat1,rat2)
        enddo
c
        return      
      END SUBROUTINE     
c
c
c
      SUBROUTINE shf_init(nterms, ntheta, nphi, theta,
     $ ctheta, whts, ynms, wsave, ylmGrid)
      implicit real *8 (a-h,o-z)
      !! input
      integer :: nterms, nphi, ntheta
      !! output
      real*8 :: theta(ntheta)
      complex *16 ylmGrid(nphi,ntheta,0:nterms,-nterms:nterms)     
      real* 8 :: ctheta(ntheta), whts(ntheta), beta(ntheta)
      real* 8 :: ynms(0:nterms,0:nterms,ntheta/2+1)
      complex *16 :: wsave(4*nphi+15)

      complex *16, ALLOCATABLE:: mpole(:,:)    
            
      allocate(mpole(0:nterms,-nterms:nterms))     

      !! prepare spherical harmonic functions for nterms
      call sphtrans_cmpl_lege_init(nterms,nphi,ntheta,ctheta,
     $ whts,ynms,wsave)

      do i = 1,ntheta
          theta(i) = dacos(ctheta(i))
      enddo     
      
      do n = 0, nterms
        do m = -n, n
          mpole(n,m) = (0d0, 0d0)
        enddo
      enddo     
            
      do n = 0, nterms
        do m = -n, n
          mpole(n,m) = (1d0, 0d0)
         ! write(*,*) "create sph grid", n, m
          call sphtrans_cmpl(nterms,mpole,nphi,ntheta,ylmGrid(:,:,n,m),
     $     ctheta,ynms,wsave)      
          mpole(n,m) = (0d0, 0d0)
        enddo
      enddo      

      END SUBROUTINE


      SUBROUTINE shf_first_dev(nterms, ntheta, nphi, 
     $ ylmGrid, dydthGrid, dydphGrid)
      implicit real *8 (a-h,o-z)
      complex *16 eiphi, conj_eiphi, tmp1, tmp2, tmp3
      !! input      
      complex *16 ylmGrid(nphi,ntheta,0:nterms,-nterms:nterms)
      !! output
      complex *16 dydthGrid(nphi,ntheta,0:nterms,-nterms:nterms)
      complex *16 dydphGrid(nphi,ntheta,0:nterms,-nterms:nterms)

!! Our definition of complex spherical harmonics is
!! Ynm(theta,phi)= sqrt(2n+1) sqrt((n-m)!/(n+m)!)  Pnm(cos theta) e^(im phi), 
!! Yn,-m(theta,phi) = sqrt(2n+1) sqrt((n-m)!/(n+m)!) Pnm(cos theta) e^(-imphi), for m >= 0.
!! Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.
      
      PI = 4.D0*DATAN(1.D0)
      phi = 2*PI/nphi
      do i=1,ntheta
        do j=1,nphi
          phij = (j-1)*phi
          eiphi = DCMPLX(dcos(phij),dsin(phij))  !! exp(i*phi)
          conj_eiphi = DCMPLX(dcos(phij),-dsin(phij))
          do l = 0,nterms
            do m = -l, l
              dydphGrid(j,i,l,m)=ylmGrid(j,i,l,m)*DCMPLX(0,m)
              tmp1 = DCMPLX(0,0)
              tmp2 = DCMPLX(0,0)
              IF (m-1 .ge. -l) then
              tmp1=dsqrt(dble((l+m)*(l-m+1)))*ylmGrid(j,i,l,m-1)*eiphi
                IF (m-1.lt.0) then
                  tmp1 = tmp1*dble((-1)**(m-1))
                ENDIF
              ENDIF
              IF (m+1 .le. l) then
           tmp2=dsqrt(dble((l+m+1)*(l-m)))*ylmGrid(j,i,l,m+1)*conj_eiphi
                IF (m+1.lt.0) then
                  tmp2 = tmp2*dble((-1)**(m+1))
                ENDIF
              ENDIF              
              dydthGrid(j,i,l,m)= -0.5d0*(tmp1 - tmp2)              
              IF (m.lt.0) then
                dydthGrid(j,i,l,m)=dydthGrid(j,i,l,m)*dble((-1)**(m))
              ENDIF
            enddo
          enddo
        enddo
      enddo     
              
      END SUBROUTINE
     
     
      SUBROUTINE upsamp_zero_pad(npterms,nterms,xmpole,ympole,zmpole,
     $ xqmpole, yqmpole, zqmpole) 
      implicit real *8 (a-h,o-z)
      integer :: npterms, nterms, m, n
      complex *16 :: xmpole(0:npterms,-npterms:npterms),
     $ ympole(0:npterms,-npterms:npterms),
     $ zmpole(0:npterms,-npterms:npterms)     
      complex *16 :: xqmpole(0:nterms,-nterms:nterms),
     $ yqmpole(0:nterms,-nterms:nterms),zqmpole(0:nterms,-nterms:nterms)
                
      !! zero padding x,y,zqmpole from pterms to nterms
      do n = 0,npterms
        do m = -n,n
          xqmpole(n,m) = xmpole(n, m)
          yqmpole(n,m) = ympole(n, m)
          zqmpole(n,m) = zmpole(n, m)                   
        enddo
      enddo      
      do n = npterms+1,nterms
        do m = -n,n
          xqmpole(n,m) = (0d0, 0d0)
          yqmpole(n,m) = (0d0, 0d0)
          zqmpole(n,m) = (0d0, 0d0)                    
        enddo
      enddo      
           
      END SUBROUTINE


      SUBROUTINE upsamp_zero_pad_real(npterms,nterms,
     $ xmpole,ympole,zmpole,xqmpole, yqmpole, zqmpole) 
      implicit real *8 (a-h,o-z)
      integer :: npterms, nterms, m, n
      complex *16 :: xmpole(0:npterms,0:npterms),
     $ ympole(0:npterms,0:npterms),
     $ zmpole(0:npterms,0:npterms)     
      complex *16 :: xqmpole(0:nterms,0:nterms),
     $ yqmpole(0:nterms,0:nterms),zqmpole(0:nterms,0:nterms)
                
      !! zero padding x,y,zqmpole from pterms to nterms
      do n = 0,npterms
        do m = 0,n
          xqmpole(n,m) = xmpole(n, m)
          yqmpole(n,m) = ympole(n, m)
          zqmpole(n,m) = zmpole(n, m)                   
        enddo
      enddo      
      do n = npterms+1,nterms
        do m = 0,n
          xqmpole(n,m) = (0d0, 0d0)
          yqmpole(n,m) = (0d0, 0d0)
          zqmpole(n,m) = (0d0, 0d0)                    
        enddo
      enddo      
           
      END SUBROUTINE

      SUBROUTINE trunc_mpole(npterms,nterms,xmpole,ympole,zmpole,
     $ xqmpole, yqmpole, zqmpole) 
      implicit real *8 (a-h,o-z)
      integer :: npterms, nterms, m, n
      complex *16 :: xmpole(0:npterms,-npterms:npterms),
     $ ympole(0:npterms,-npterms:npterms),
     $ zmpole(0:npterms,-npterms:npterms)     
      complex *16 :: xqmpole(0:nterms,-nterms:nterms),
     $ yqmpole(0:nterms,-nterms:nterms),zqmpole(0:nterms,-nterms:nterms)
                
      !! zero padding x,y,zqmpole from pterms to nterms
      do n = 0,npterms
        do m = -n,n
          xmpole(n, m) = xqmpole(n,m) 
          ympole(n, m) = yqmpole(n,m) 
          zmpole(n, m) = zqmpole(n,m)             
        enddo
      enddo      
      END SUBROUTINE
