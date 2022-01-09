!! input:
!!   nterms, ntheta, nphi,
!!   whts, theta,
!!   Geom_W, xgrid, ygrid, zgrid
!! output:
!!   surf_area, tot_vol
      SUBROUTINE eval_area_vol(ntheta, nphi, theta, whts,
     $ xgrid, ygrid, zgrid, Geom_W, norm, surf_area, tot_vol)

      implicit real *8 (a-h,o-z)
        integer :: ntheta, nphi
        real* 8 :: theta(ntheta), whts(ntheta)
        complex*16 :: xgrid(nphi,ntheta), ygrid(nphi,ntheta),
     $ zgrid(nphi,ntheta)
        real* 8 :: Geom_W(nphi,ntheta), norm(nphi, ntheta, 3)  
        real* 8 :: surf_area, tot_vol
        real* 8 :: radius
        
        PI = 4.D0*DATAN(1.D0)

        tmp2 = 0d0
        do i = 1, ntheta
            tmp1 = 0d0
            do j = 1, nphi
                !! "/dsin(theta)" because the integrated variable is
                !! not dtheta, but dx, where x = cos(theta) from -1 to 1
                tmp1 = tmp1 + Geom_W(j,i)/dsin(theta(i))
            enddo         
            tmp2 = tmp2 + tmp1 * whts(i)  
        enddo
        surf_area = tmp2/dble(nphi)*2*pi
        write(100,*) "surf_area = ", surf_area

        tmp2 = 0d0
        do i = 1, ntheta
            tmp1 = 0d0
            do j = 1, nphi
                radius = dsqrt(dble(xgrid(j,i))**2 + dble(ygrid(j,i))**2
     $        + dble(zgrid(j,i))**2)
                !write(*,*) "radius = ", radius 
                tmp = norm(j,i,1)*dble(xgrid(j,i))
     $ + norm(j,i,2)*dble(ygrid(j,i)) + norm(j,i,3)*dble(zgrid(j,i))
                tmp1 = tmp1 + Geom_W(j,i) * tmp / 3d0 /dsin(theta(i))
            enddo
            tmp2 = tmp2 + tmp1 * whts(i) 
        enddo
        tot_vol = tmp2/dble(nphi)*2*pi
        write(100,*) "tot_vol = ", tot_vol
        
      END SUBROUTINE
c
c
c
      SUBROUTINE eval_area_vol_real(ntheta, nphi, theta, whts,
     $ xgrid, ygrid, zgrid, Geom_W, norm, surf_area, tot_vol)
   
         implicit real *8 (a-h,o-z)
           integer :: ntheta, nphi
           real* 8 :: theta(ntheta), whts(ntheta)
           real* 8 :: xgrid(nphi,ntheta), ygrid(nphi,ntheta),
     $ zgrid(nphi,ntheta)
           real* 8 :: Geom_W(nphi,ntheta), norm(nphi, ntheta, 3)  
           real* 8 :: surf_area, tot_vol
           real* 8 :: radius
           
           PI = 4.D0*DATAN(1.D0)
   
           tmp2 = 0d0
           do i = 1, ntheta
               tmp1 = 0d0
               do j = 1, nphi
                   !! "/dsin(theta)" because the integrated variable is
                   !! not dtheta, but dx, where x = cos(theta) from -1 to 1
                   tmp1 = tmp1 + Geom_W(j,i)/dsin(theta(i))
               enddo         
               tmp2 = tmp2 + tmp1 * whts(i)  
           enddo
           surf_area = tmp2/dble(nphi)*2*pi
   
           tmp2 = 0d0
           do i = 1, ntheta
               tmp1 = 0d0
               do j = 1, nphi
                   radius = dsqrt(xgrid(j,i)**2 + ygrid(j,i)**2
     $        + zgrid(j,i)**2)
                   !write(*,*) "radius = ", radius 
                   tmp = norm(j,i,1) * xgrid(j,i)
     $  + norm(j,i,2) * ygrid(j,i) + norm(j,i,3) * zgrid(j,i)
                   tmp1 = tmp1 + Geom_W(j,i) * tmp / 3d0 /dsin(theta(i))
               enddo
               tmp2 = tmp2 + tmp1 * whts(i) 
           enddo
           tot_vol = tmp2 / dble(nphi) * 2 * pi
           write(100,*) "surf_area = ", surf_area, "tot_vol = ", tot_vol
           
         END SUBROUTINE
c
c
c
      SUBROUTINE eval_area_vol_cenMass_real(ntheta, nphi, whts,
     $ xgrid, ygrid, zgrid, wsgrid, norm, surf_area, tot_vol, cen_mass)
   
         implicit real *8 (a-h,o-z)
           integer :: ntheta, nphi
           real* 8 :: whts(ntheta)
           real* 8 :: xgrid(nphi,ntheta), ygrid(nphi,ntheta),
     $ zgrid(nphi,ntheta)
           real* 8 :: wsgrid(nphi,ntheta), norm(nphi, ntheta, 3) 
           real *8 :: cen_mass(3) 
           real* 8 :: surf_area, tot_vol
           real* 8 :: radius
           
           PI = 4.D0*DATAN(1.D0)
   
           tmp2 = 0d0
           do i = 1, ntheta
               tmp1 = 0d0
               do j = 1, nphi
                   !! "/dsin(theta)" because the integrated variable is
                   !! not dtheta, but dx, where x = cos(theta) from -1 to 1
                   tmp1 = tmp1 + wsgrid(j,i)
               enddo         
               tmp2 = tmp2 + tmp1 * whts(i)  
           enddo
           surf_area = tmp2/dble(nphi)*2*pi
   
           tmp2 = 0d0
           do i = 1, ntheta
               tmp1 = 0d0
               do j = 1, nphi
                   !write(*,*) "radius = ", radius 
                   tmp = norm(j,i,1) * xgrid(j,i)
     $  + norm(j,i,2) * ygrid(j,i) + norm(j,i,3) * zgrid(j,i)
                   tmp1 = tmp1 + wsgrid(j,i) * tmp / 3d0 
               enddo
               tmp2 = tmp2 + tmp1 * whts(i) 
           enddo
           tot_vol = tmp2 / dble(nphi) * 2 * pi

           tmp2x = 0d0
           tmp2y = 0d0
           tmp2z = 0d0
           do i = 1, ntheta
               tmp1x = 0d0
               tmp1y = 0d0
               tmp1z = 0d0
               do j = 1, nphi
                   !write(*,*) "radius = ", radius 
                   tmp = norm(j,i,1) * xgrid(j,i)
     $  + norm(j,i,2) * ygrid(j,i) + norm(j,i,3) * zgrid(j,i)
                   tmp1x = tmp1x + xgrid(j,i) * wsgrid(j,i) * tmp / 4d0 
                   tmp1y = tmp1y + ygrid(j,i) * wsgrid(j,i) * tmp / 4d0 
                   tmp1z = tmp1z + zgrid(j,i) * wsgrid(j,i) * tmp / 4d0 
               enddo
               tmp2x = tmp2x + tmp1x * whts(i) 
               tmp2y = tmp2y + tmp1y * whts(i) 
               tmp2z = tmp2z + tmp1z * whts(i) 
           enddo
           cen_mass(1) = tmp2x / dble(nphi) * 2 * pi / tot_vol 
           cen_mass(2) = tmp2y / dble(nphi) * 2 * pi / tot_vol 
           cen_mass(3) = tmp2z / dble(nphi) * 2 * pi / tot_vol 

         END SUBROUTINE
