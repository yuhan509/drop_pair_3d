      SUBROUTINE close_dist_flag(iflag, x0, y0, z0,
     $ xgrid, ygrid, zgrid, xmpole, ympole, zmpole,
     $ dist_cl, theta_cl, phi_cl, sphf_s, xs, ys, zs)

      use mod_orgHelper       
      implicit real *8 (a-h,o-z)

      real *8 :: x0, y0, z0, dist_cl
      integer :: iflag   !! disable search when iflag > 1, return closest grid point

      complex *16 :: xmpole(0:nterms,0:nterms),
     $ ympole(0:nterms,0:nterms), zmpole(0:nterms,0:nterms)
      real *8 :: xgrid(nphi,ntheta),ygrid(nphi,ntheta),
     $ zgrid(nphi,ntheta)      
      complex *16 :: sphf_s(0:nterms, -nterms:nterms)
      real *8 :: xs, ys, zs, theta_cl, phi_cl


      !! find closest grid point (j_cl, i_cl)
      dist_min = 1d11
      do i = 1, ntheta
         do j = 1, nphi
            dist = (x0 - xgrid(j,i))**2 
     $           + (y0 - ygrid(j,i))**2
     $           + (z0 - zgrid(j,i))**2
            if (dist_min .gt. dist) then
              dist_min = dist
              j_cl = j
              i_cl = i
            endif
         enddo
      enddo
      dist_min = sqrt(dist_min)
      
      if (iflag .gt. 1) then
        dist_cl = dist_min
        return
      endif

      pi = 4*datan(1d0)
   !!!! ****** debugging disable search_closest_phi_theta
       if (dist_min .lt. pi/(nterms+1)*2) then 
    !  if (dist_min .lt. 6 * pi/(nterms+1)) then  !!! at long distance no need for dist_cl evaluation
   !!!!

   !   write(910, *) "closest grid point ", "j =", j_cl, "i =", i_cl, 
   !  $ "phi =",phi_cl,"theta =", theta_cl
      !! use x,y,zgrid(j,i) as the starting point for Newton's Method
      !! find the closest point (phi_cl, theta_cl) on the surface
        
      call search_closest_phi_theta(x0, y0, z0, j_cl, i_cl, 
     $ xmpole, ympole, zmpole, dist_cl, theta_cl, phi_cl,
     $ sphf_s, xs, ys, zs)

      if (dist_cl .lt. 0d0) then   !! negative determinant of hessian
          dist_cl = dist_min
      endif
   
      else
       dist_cl = dist_min
      endif
   
      END SUBROUTINE
c
c
c
      SUBROUTINE search_closest_phi_theta(x0, y0, z0, j_cl, i_cl, 
     $ xmpole, ympole, zmpole, dist_cl, theta_cl, phi_cl,
     $ sphf, xs, ys, zs)

!      use mod_upsamp
      use mod_orgHelper
      implicit real *8 (a-h,o-z)      

      real *8 :: x0, y0, z0, theta_cl, phi_cl, dist_cl
      real *8 :: xs, ys, zs
      complex *16 :: xmpole(0:nterms,0:nterms),
     $ ympole(0:nterms,0:nterms), zmpole(0:nterms,0:nterms)   
      complex *16 :: sphf(0:nterms, -nterms:nterms)

      complex *16, ALLOCATABLE, dimension(:,:) :: 
     $ dydph, dydth, d2ydth2, d2ydph2, d2ydthdph  
      complex *16 :: ctmp
      allocate(dydth(0:nterms,-nterms:nterms))
      allocate(dydph(0:nterms,-nterms:nterms))
      allocate(d2ydth2(0:nterms,-nterms:nterms))
      allocate(d2ydthdph(0:nterms,-nterms:nterms))
      allocate(d2ydph2(0:nterms,-nterms:nterms))
         
      pi = 4d0 * datan(1d0)
      phi_cl = 2*pi / nphi * (j_cl - 1)
      theta_cl = theta(i_cl) 
      
    ! write(920, *) "start: phi_cl = ", phi_cl, "theta_cl", theta_cl
      tol = 1d-12

      iloop = 0
 20   continue

       call sphf_and_derivs_atapoint(nterms, theta_cl, phi_cl,
     $ sphf, dydph, dydth, d2ydth2, d2ydph2, d2ydthdph) 

!!!   test derivatives
!      iup = 1
!      iterms = nterms
!      j = j_cl
!      i = i_cl
!      do n = 0, iterms
!         do m = -n, n
!      !     ctmp = up_sphfgrid(iup,j,i,n,m) - sphf(n,m) 
!      !     ctmp =up_dydph(iup,j,i,n,m) - dydph(n,m)
!      !     ctmp =up_dydth(iup,j,i,n,m) - dydth(n,m)
!           ctmp =up_d2ydth2(iup,j,i,n,m) - d2ydth2(n,m)
!       !    ctmp =up_d2ydph2(iup,j,i,n,m) - d2ydph2(n,m)
!       !    ctmp =up_d2ydthdph(iup,j,i,n,m) - d2ydthdph(n,m)
!           if (dble(ctmp)**2 + imag(ctmp)**2 .gt. 1d-20) then
!            write(910,*)  n, m, ctmp
!           endif
!         enddo
!      enddo

        tmpx = 0d0
        tmpx1 = 0d0
        tmpx2 = 0d0
        tmpx11 = 0d0
        tmpx12 = 0d0
        tmpx22 = 0d0
        tmpy = 0d0
        tmpy1 = 0d0
        tmpy2 = 0d0
        tmpy11 = 0d0
        tmpy12 = 0d0
        tmpy22 = 0d0
        tmpz = 0d0
        tmpz1 = 0d0
        tmpz2 = 0d0
        tmpz11 = 0d0
        tmpz12 = 0d0
        tmpz22 = 0d0
        do n = 0, nterms
            m = 0
            tmpx = tmpx + xmpole(n,m) * sphf(n,m)
            tmpx1 = tmpx1 + xmpole(n,m) * dydph(n,m)
            tmpx2 = tmpx2 + xmpole(n,m) * dydth(n,m)
            tmpx11 = tmpx11 + xmpole(n,m) * d2ydph2(n,m)
            tmpx12 = tmpx12 + xmpole(n,m) * d2ydthdph(n,m)
            tmpx22 = tmpx22 + xmpole(n,m) * d2ydth2(n,m)
            tmpy = tmpy + ympole(n,m) * sphf(n,m)
            tmpy1 = tmpy1 + ympole(n,m) * dydph(n,m)
            tmpy2 = tmpy2 + ympole(n,m) * dydth(n,m)
            tmpy11 = tmpy11 + ympole(n,m) * d2ydph2(n,m)
            tmpy12 = tmpy12 + ympole(n,m) * d2ydthdph(n,m)
            tmpy22 = tmpy22 + ympole(n,m) * d2ydth2(n,m)
            tmpz = tmpz + zmpole(n,m) * sphf(n,m)
            tmpz1 = tmpz1 + zmpole(n,m) * dydph(n,m)
            tmpz2 = tmpz2 + zmpole(n,m) * dydth(n,m)
            tmpz11 = tmpz11 + zmpole(n,m) * d2ydph2(n,m)
            tmpz12 = tmpz12 + zmpole(n,m) * d2ydthdph(n,m)
            tmpz22 = tmpz22 + zmpole(n,m) * d2ydth2(n,m)
            do m = 1, n
               tmpx = tmpx + 2*dble(xmpole(n,m) * sphf(n,m))
               tmpx1 = tmpx1 + 2*dble(xmpole(n,m) * dydph(n,m))
               tmpx2 = tmpx2 + 2*dble(xmpole(n,m) * dydth(n,m))
               tmpx11 = tmpx11 + 2*dble(xmpole(n,m) * d2ydph2(n,m))
               tmpx12 = tmpx12 + 2*dble(xmpole(n,m) * d2ydthdph(n,m))
               tmpx22 = tmpx22 + 2*dble(xmpole(n,m) * d2ydth2(n,m))
               tmpy = tmpy + 2*dble(ympole(n,m) * sphf(n,m))
               tmpy1 = tmpy1 + 2*dble(ympole(n,m) * dydph(n,m))
               tmpy2 = tmpy2 + 2*dble(ympole(n,m) * dydth(n,m))
               tmpy11 = tmpy11 + 2*dble(ympole(n,m) * d2ydph2(n,m))
               tmpy12 = tmpy12 + 2*dble(ympole(n,m) * d2ydthdph(n,m))
               tmpy22 = tmpy22 + 2*dble(ympole(n,m) * d2ydth2(n,m))
               tmpz = tmpz + 2*dble(zmpole(n,m) * sphf(n,m))
               tmpz1 = tmpz1 + 2*dble(zmpole(n,m) * dydph(n,m))
               tmpz2 = tmpz2 + 2*dble(zmpole(n,m) * dydth(n,m))
               tmpz11 = tmpz11 + 2*dble(zmpole(n,m) * d2ydph2(n,m))
               tmpz12 = tmpz12 + 2*dble(zmpole(n,m) * d2ydthdph(n,m))
               tmpz22 = tmpz22 + 2*dble(zmpole(n,m) * d2ydth2(n,m))
            enddo
        enddo

      !! Gradient evaluation
        g1 = (tmpx - x0) * tmpx1 + 
     $     (tmpy - y0) * tmpy1 + (tmpz - z0) * tmpz1
        g2 = (tmpx - x0) * tmpx2 + 
     $     (tmpy - y0) * tmpy2 + (tmpz - z0) * tmpz2
        g1 = 2*g1
        g2 = 2*g2
      !! Hessian evaluation
        h11 = tmpx1**2 + tmpy1**2 + tmpz1**2 +
     $   (tmpx - x0)*tmpx11 + (tmpy - y0)*tmpy11 + (tmpz - z0)*tmpz11
        h11 = 2 * h11
        h22 = tmpx2**2 + tmpy2**2 + tmpz2**2 +
     $   (tmpx - x0)*tmpx22 + (tmpy - y0)*tmpy22 + (tmpz - z0)*tmpz22
        h22 = 2 * h22
        h12 = tmpx1*tmpx2 + tmpy1*tmpy2 + tmpz1*tmpz2 +
     $   (tmpx - x0)*tmpx12 + (tmpy - y0)*tmpy12 + (tmpz - z0)*tmpz12
        h12 = 2 * h12      

        if (h12*h12 - h11*h22 .ge. 0d0) then 
!            tmpmax = max(h11 + abs(h12), h22 + abs(h12))
!            tmpmin = min(h11 - abs(h12), h22 - abs(h12))
!            tmp = max(tmpmax , abs(tmpmin))
!            h11 = h11 + tmp
!            h22 = h22 + tmp
            write(*,*) "WARNING: h12*h12 - h11*h22 .ge. 0d0"
            write(*,*) "j_cl,i_cl,",j_cl,i_cl,"iloop=",iloop
            dist_cl = -1d0
        else 

        !! invert Hessian (or modified Hessian)
        tmp = h11*h22 - h12*h12 
        h11i =  h22 / tmp
        h12i =  - h12 / tmp
        h22i =  h11 / tmp

        !! Newton iteration method
        !write(920, *) "h11i",h11i,"h12i",h12i,"h22i",h22i,"g1",g1,"g2",g2
        phi_cl = phi_cl - (h11i * g1 + h12i * g2)
        theta_cl = theta_cl - (h12i * g1 + h22i * g2)
    !   write(920, *) "phi_cl = ", phi_cl, "theta_cl", theta_cl

        if (((abs(h11i * g1 + h12i * g2) .gt. tol) .or. 
     $ (abs(h12i * g1 + h22i * g2) .gt. tol)) .and. iloop .lt. 25) then
    !   write(920,*) "step_phi =", h11i * g1 + h12i * g2
    !   write(920,*) "step_theta =", h12i * g1 + h22i * g2
    !   write(920,*) "continue", iloop
          iloop = iloop + 1
          goto 20
        endif

cccccccccc eval dist cccccccc
        call sphf_and_derivs_atapoint(nterms, theta_cl, phi_cl,
     $ sphf, dydph, dydth, d2ydth2, d2ydph2, d2ydthdph) 

        tmpx = 0d0
        tmpy = 0d0
        tmpz = 0d0
        do n = 0, nterms
            m = 0
            tmpx = tmpx + xmpole(n,m) * sphf(n,m)
            tmpy = tmpy + ympole(n,m) * sphf(n,m)
            tmpz = tmpz + zmpole(n,m) * sphf(n,m)
            do m = 1, n
               tmpx = tmpx + 2*dble(xmpole(n,m) * sphf(n,m))
               tmpy = tmpy + 2*dble(ympole(n,m) * sphf(n,m))
               tmpz = tmpz + 2*dble(zmpole(n,m) * sphf(n,m))
            enddo
        enddo
        xs = tmpx
        ys = tmpy
        zs = tmpz

        dist = (x0 - xs)**2 + (y0 - ys)**2 + (z0 - zs)**2
        dist_cl = dsqrt(dist)
    !    write(920, *) "dist_cl = ", dist_cl
cccccccccccccccccccccc
      endif
      END SUBROUTINE
c
c
c
      SUBROUTINE sphf_and_derivs_atapoint(nterms, theta_s, phi_s,
     $ sphf, dydph, dydth, d2ydth2, d2ydph2, d2ydthdph) 

       implicit real *8 (a-h,o-z)

       integer :: nterms
       real *8 :: theta_s, phi_s
       complex *16 eiphi, conj_eiphi, tmp1, tmp2, tmp3      
       complex *16 :: sphf(0:nterms,-nterms:nterms), 
     $ dydph(0:nterms,-nterms:nterms), dydth(0:nterms,-nterms:nterms),
     $ d2ydth2(0:nterms,-nterms:nterms), 
     $ d2ydph2(0:nterms,-nterms:nterms), 
     $ d2ydthdph(0:nterms,-nterms:nterms)

       real *8, allocatable :: rat1(:,:),rat2(:,:),ynm(:,:)

       allocate(rat1(nterms+1,nterms+1))
       allocate(rat2(nterms+1,nterms+1))
       allocate(ynm(0:nterms,0:nterms))

c     subroutine ylgndrini(nmax, rat1, rat2)
c     Precompute the recurrence coefficients for the fast
c     evaluation of normalized Legendre functions and their derivatives
c    
c     Parameters:
c       nmax                      must be non-negative
c       rat1(0:nmax,0:nmax)       recurrence coefficient
c       rat2(0:nmax,0:nmax)       recurrence coefficient
       call ylgndrini(nterms,rat1,rat2) 

c     subroutine ylgndrf(nmax, x, y, rat1, rat2)       
c     Evaluate normalized Legendre functions
c
c      Ynm(x) = sqrt(2n+1)  sqrt( (n-m)!/ (n+m)! ) Pnm(x)
c
c     for n = 0, 1, 2,..., nmax
c     and  m = 0, 1,..., n.
c
c     Parameters:
c     nmax                  must be non-negative
c     x                     -1 <= x <= 1
c     y(0:nmax,0:nmax)      resulting function values
c
c     Upon return, y(n,m) will contain the function value Ynm(x) for 0
c     <= n <= nmax and 0 <= m <= n.  Other elements of y will contain
c     undefined values.             

    !     write(910, *) phi_s, theta_s
          call ylgndrf(nterms,dcos(theta_s),ynm,rat1,rat2)
          do n = 0, nterms
            do m = -n, n
              sphf(n,m) = ynm(n,abs(m))
     $         *DCMPLX(dcos(m*phi_s),dsin(m*phi_s))
            enddo
          enddo   

          phij = phi_s
          eiphi = DCMPLX(dcos(phij),dsin(phij))  !! exp(i*phi)
          conj_eiphi = DCMPLX(dcos(phij),-dsin(phij))
          do l = 0,nterms
            do m = -l, l
              dydph(l,m)=sphf(l,m)*DCMPLX(0,m)
              d2ydph2(l,m)=sphf(l,m)*dble(-m*m)
              tmp1 = DCMPLX(0,0)
              tmp2 = DCMPLX(0,0)
              IF (m-1 .ge. -l) then
              tmp1=dsqrt(dble((l+m)*(l-m+1)))*sphf(l,m-1)*eiphi
                IF (m-1.lt.0) then
                  tmp1 = tmp1*dble((-1)**(m-1))
                ENDIF
              ENDIF
              IF (m+1 .le. l) then
           tmp2=dsqrt(dble((l+m+1)*(l-m)))*sphf(l,m+1)*conj_eiphi
                IF (m+1.lt.0) then
                  tmp2 = tmp2*dble((-1)**(m+1))
                ENDIF
              ENDIF
              
              dydth(l,m)= -0.5d0*(tmp1 - tmp2)              
              IF (m.lt.0) then
                dydth(l,m)=dydth(l,m)*dble((-1)**(m))
              ENDIF
              
              d2ydthdph(l,m) = dydth(l,m)*DCMPLX(0,m)
            enddo
          enddo

          phij = 2*phi_s   !!! it's actually 2 * phi_s
          eiphi = DCMPLX(dcos(phij),dsin(phij))
          conj_eiphi = DCMPLX(dcos(phij),-dsin(phij))
          do l = 0,nterms
            do m = -l, l              
              tmp1 = DCMPLX(0,0)
              tmp2 = DCMPLX(0,0)
              tmp3 = dble(2*(l**2-m**2+l))*sphf(l,m)
              IF (m .lt. 0) then
                tmp3 = tmp3 * dble((-1)**m)
              ENDIF
              IF (m-2 .ge. -l) then
              tmp1=dsqrt(dble((l+m)*(l-m+1)*(l+m-1)*(l-m+2)))
     $ *sphf(l,m-2)*eiphi
                 IF (m-2.lt.0) then
                  tmp1 = tmp1*dble((-1)**(m-2))
                 ENDIF
              ENDIF
              IF (m+2 .le. l) then
               tmp2=dsqrt(dble((l+m+1)*(l-m))*(l+m+2)*(l-m-1))
     $ *sphf(l,m+2)*conj_eiphi
                IF (m+2.lt.0) then
                  tmp2 = tmp2*dble((-1)**(m+2))
                ENDIF
              ENDIF
              d2ydth2(l,m)= 0.25d0*(tmp1 + tmp2 - tmp3)              
              IF (m.lt.0) then
                d2ydth2(l,m) = d2ydth2(l,m)*dble((-1)**(m))
              ENDIF           
            enddo
          enddo

      END SUBROUTINE


      SUBROUTINE close_dist_test(x0, y0, z0,
     $ xgrid, ygrid, zgrid, xmpole, ympole, zmpole,
     $ theta_cl, phi_cl, dist_cl)

      use mod_orgHelper       
      implicit real *8 (a-h,o-z)

      real *8 :: x0, y0, z0, dist_cl

      complex *16 :: xmpole(0:nterms,0:nterms),
     $ ympole(0:nterms,0:nterms), zmpole(0:nterms,0:nterms)
      real *8 :: xgrid(nphi,ntheta),ygrid(nphi,ntheta),
     $ zgrid(nphi,ntheta)      
      complex *16 :: sphf_s(0:nterms, -nterms:nterms)
      real *8 :: xs, ys, zs


      !! find closest grid point (j_cl, i_cl)
      dist_min = 1d11
      do i = 1, ntheta
         do j = 1, nphi
            dist = (x0 - xgrid(j,i))**2 
     $           + (y0 - ygrid(j,i))**2
     $           + (z0 - zgrid(j,i))**2
            if (dist_min .gt. dist) then
              dist_min = dist
              j_cl = j
              i_cl = i
            endif
         enddo
      enddo


   !   write(910, *) "closest grid point ", "j =", j_cl, "i =", i_cl, 
   !  $ "phi =",phi_cl,"theta =", theta_cl
      !! use x,y,zgrid(j,i) as the starting point for Newton's Method
      !! find the closest point (phi_cl, theta_cl) on the surface
      call search_closest_phi_theta(x0, y0, z0, j_cl, i_cl, 
     $ xmpole, ympole, zmpole, dist_cl,
     $ sphf_s, xs, ys, zs)
   
      END SUBROUTINE
