      SUBROUTINE shf_diffs(ylmGrid,nterms,nphi,ntheta,dydthGrid, 
     $ dydphGrid,d2ydth2Grid,d2ydthdphGrid,d2ydph2Grid)
      implicit real *8 (a-h,o-z)
      complex *16 eiphi, conj_eiphi, tmp1, tmp2, tmp3      
      complex *16 ylmGrid(nphi,ntheta,0:nterms,-nterms:nterms)
      complex *16 dydthGrid(nphi,ntheta,0:nterms,-nterms:nterms)
      complex *16 dydphGrid(nphi,ntheta,0:nterms,-nterms:nterms)
      complex *16 d2ydth2Grid(nphi,ntheta,0:nterms,-nterms:nterms)
      complex *16 d2ydthdphGrid(nphi,ntheta,0:nterms,-nterms:nterms)            
      complex *16 d2ydph2Grid(nphi,ntheta,0:nterms,-nterms:nterms)

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
              d2ydph2Grid(j,i,l,m)=ylmGrid(j,i,l,m)*dble(-m*m)
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
              
              d2ydthdphGrid(j,i,l,m) = dydthGrid(j,i,l,m)*DCMPLX(0,m)
            enddo
          enddo
        enddo
      enddo

      
      do i=1,ntheta
        do j=1,nphi
          phij = 2*(j-1)*phi   !!! it's actually 2 * phij
          eiphi = DCMPLX(dcos(phij),dsin(phij))
          conj_eiphi = DCMPLX(dcos(phij),-dsin(phij))
          do l = 0,nterms
            do m = -l, l
            
!              IF ((m+2 .le. l) .and. (m-2 .ge. -l)) then
!                d2ydth2Grid(j,i,l,m) = 0.25d0*
!     $ (dsqrt(dble((l+m)*(l-m+1)*(l+m-1)*(l-m+2)))
!     $ *ylmGrid(j,i,l,m-2)*eiphi
!     $ - dble(2*(n**2-m**2+n))*ylmGrid(j,i,l,m)
!     $ + dsqrt(dble((l+m+1)*(l-m))*(l+m+2)*(l-m-1))
!     $ *ylmGrid(j,i,l,m+2)*conj_eiphi)
!              ELSE IF (m+2 .gt. l) then 
!                d2ydth2Grid(j,i,l,m) = 0.25d0*
!     $ (dsqrt(dble((l+m)*(l-m+1)*(l+m-1)*(l-m+2)))
!     $ *ylmGrid(j,i,l,m-2)*eiphi
!     $ - dble(2*(n**2-m**2+n))*ylmGrid(j,i,l,m))
!              ELSE IF (m-2 .lt. -l) then
!                d2ydth2Grid(j,i,l,m) = 0.25d0*
!     $ (-dble(2*(n**2-m**2+n))*ylmGrid(j,i,l,m)
!     $ + dsqrt(dble((l+m+1)*(l-m))*(l+m+2)*(l-m-1))
!     $ *ylmGrid(j,i,l,m+2)*conj_eiphi)
!              ENDIF
              
              tmp1 = DCMPLX(0,0)
              tmp2 = DCMPLX(0,0)
              tmp3 = dble(2*(l**2-m**2+l))*ylmGrid(j,i,l,m)
              IF (m .lt. 0) then
                tmp3 = tmp3 * dble((-1)**m)
              ENDIF
!              IF (l.eq.1) then
!                write(27,*) j,i,l,m,'tmp3',tmp3,'ylm',ylmGrid(j,i,l,m)
!              ENDIF  
              IF (m-2 .ge. -l) then
              tmp1=dsqrt(dble((l+m)*(l-m+1)*(l+m-1)*(l-m+2)))
     $ *ylmGrid(j,i,l,m-2)*eiphi
                 IF (m-2.lt.0) then
                  tmp1 = tmp1*dble((-1)**(m-2))
                 ENDIF
              ENDIF
              IF (m+2 .le. l) then
               tmp2=dsqrt(dble((l+m+1)*(l-m))*(l+m+2)*(l-m-1))
     $ *ylmGrid(j,i,l,m+2)*conj_eiphi
                IF (m+2.lt.0) then
                  tmp2 = tmp2*dble((-1)**(m+2))
                ENDIF
              ENDIF
              d2ydth2Grid(j,i,l,m)= 0.25d0*(tmp1 + tmp2 - tmp3)              
              IF (m.lt.0) then
                d2ydth2Grid(j,i,l,m)=
     $ d2ydth2Grid(j,i,l,m)*dble((-1)**(m))
              ENDIF 
              IF ((j.eq.6) .and. (i.eq.1)) then
!              write(27,*) j,i,l,m,'ylm',ylmGrid(j,i,l,m),
!     $    'tmp1',tmp1,'tmp2',tmp2,'tmp3',tmp3,
!     $      d2ydth2Grid(j,i,l,m)    
              ENDIF                
            enddo
          enddo
        enddo
      enddo      
            
      END SUBROUTINE
