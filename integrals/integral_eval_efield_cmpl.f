       SUBROUTINE nearSing_integral_efield_cmpl(
     $ iterms, itheta, iphi,
     $ xigrid, yigrid, zigrid, wsigrid, sphfigrid, whts_i,
     $ gxmpole, gympole, gzmpole, dist_cl, sphf,  !! interpolate singluar integral
     $ xs, ys, zs, x0, y0, z0, norm0, cres)
      
      use mod_orgHelper     !! need to change name: ctheta... 

      implicit real *8 (a-h,o-z)
      integer :: iterms, itheta, iphi
      real *8 :: x0, y0, z0, xs, ys, zs
      real *8 :: norm0(3)
      complex *16 :: cres
      !! needed in integrand
      real *8 :: xigrid(iphi,itheta), yigrid(iphi,itheta),
     $ zigrid(iphi,itheta), wsigrid(iphi,itheta)
      real *8 :: whts_i(itheta)
      complex *16 :: sphfigrid(iphi,itheta)
      !! needed in singular integral interpolation
      complex *16 :: gxmpole(0:nterms, -nterms:nterms), 
     $ gympole(0:nterms, -nterms:nterms), 
     $ gzmpole(0:nterms, -nterms:nterms), 
     $ sphf(0:nterms, -nterms:nterms)

      real *8 :: dist_cl
      
      complex *16 :: ctmp
      real *8, ALLOCATABLE :: xarray(:)
      complex * 16, ALLOCATABLE :: farray(:)

      PI = 4.D0*DATAN(1.D0)
      h = PI/(nterms + 1)/20
      if (dist_cl .ge. h) then

        call integral_eval_efield_cmpl(iterms,itheta,iphi,
     $ xigrid, yigrid, zigrid, wsigrid, sphfigrid, whts_i,
     $ x0, y0, z0, norm0, cres)
           
       else

      allocate(xarray(nintpol),farray(nintpol))      
      !! interpolate singluar integral from expansion
      ctmp1 = (0, 0)
      ctmp2 = (0, 0)
      ctmp3 = (0, 0)
      do n = 0, nterms
        do m = -n, n
          ctmp1 = ctmp1 + gxmpole(n, m) * sphf(n, m)
          ctmp2 = ctmp2 + gympole(n, m) * sphf(n, m)
          ctmp3 = ctmp3 + gzmpole(n, m) * sphf(n, m)
        enddo
      enddo

      xarray(1) = 0d0 !! point on the surface
      farray(1) = ctmp1 * norm0(1) + ctmp2 * norm0(2) + ctmp3 * norm0(3)
!      write(903,*) 1, "d = ", xarray(1), "res = ", ctmp
      
      ratio_x = (x0 - xs) / dist_cl
      ratio_y = (y0 - ys) / dist_cl
      ratio_z = (z0 - zs) / dist_cl

!      write(903,*) "h=", h
      do i = 2, nintpol 
        d = dsqrt(dble(i-1)) * hintpol !! modified  
        xarray(i) = d + xarray(1)
        x = xs + ratio_x * d
        y = ys + ratio_y * d
        z = zs + ratio_z * d
        call integral_eval_efield_cmpl(iterms,itheta,iphi,
     $ xigrid, yigrid, zigrid, wsigrid, sphfigrid, whts_i,
     $ x, y, z, norm0, cres)  
        farray(i) = cres
 !       write(903,*) i, "d = ", d, "res = ", cres
      enddo
      
      !! target point
      xi = dist_cl

      call lagrange_value_1d_cmpl(nintpol, xarray, farray, 1, xi, cres)
!      write(903,*) "yi =", cres

       endif
       return
       END SUBROUTINE
c
c
c             
      SUBROUTINE integral_eval_efield_cmpl(iterms,itheta,iphi,
     $ xigrid, yigrid, zigrid, wsigrid, sphfgrid, whts_i,
     $ x0, y0, z0, norm0, cres)  

      implicit real *8 (a-h,o-z)
      integer :: iterms,itheta,iphi
      real *8 :: x0, y0, z0
      real *8 :: norm0(3)
      complex *16 :: cres
      real *8 :: xigrid(iphi,itheta), yigrid(iphi,itheta),
     $ zigrid(iphi,itheta), wsigrid(iphi,itheta)
      real *8 :: whts_i(itheta)
      complex *16 :: sphfgrid(iphi,itheta)
    
      complex *16 :: ctmp, cquadr

      PI = 4.D0*DATAN(1.D0)

cccccc regular gauss-legendre quadrature  ccc

        cquadr = (0d0,0d0)
        do i = 1,itheta
          ctmp = (0d0,0d0)
          do j = 1, iphi
            dist = (x0-xigrid(j,i))**2 
     $             + (y0-yigrid(j,i))**2
     $             + (z0-zigrid(j,i))**2
            dist = dsqrt(dist)
            tmp_inner = norm0(1) * (x0-xigrid(j,i))
     $          + norm0(2) * (y0-yigrid(j,i)) 
     $          + norm0(3) * (z0-zigrid(j,i)) 
            ctmp = ctmp + tmp_inner/(dist**3)*sphfgrid(j,i)
     $ * wsigrid(j,i)
  !         write(931,*) j,i, "dist,tmp_inner,ctmp"
  !         write(931,*) dist,tmp_inner,ctmp
          enddo
          cquadr = ctmp * whts_i(i) + cquadr
        enddo    
        cres = cquadr/(-2*iphi)
              
      end SUBROUTINE       
c
c
c
      SUBROUTINE integral_eval_efield_cmpl3(iterms,itheta,iphi,
     $ xigrid, yigrid, zigrid, wsigrid, sphfgrid, whts_i,
     $ x0, y0, z0, cquadr1, cquadr2, cquadr3)  

      implicit real *8 (a-h,o-z)
      integer :: iterms,itheta,iphi
      real *8 :: x0, y0, z0
      real *8 :: xigrid(iphi,itheta), yigrid(iphi,itheta),
     $ zigrid(iphi,itheta), wsigrid(iphi,itheta)
      real *8 :: whts_i(itheta)
      complex *16 :: sphfgrid(iphi,itheta)    
      complex *16 :: ctmp, ctmp1, ctmp2, ctmp3
      complex *16 :: cquadr1, cquadr2, cquadr3

      PI = 4.D0*DATAN(1.D0)

cccccc regular gauss-legendre quadrature  ccc

        cquadr1 = (0d0,0d0)
        cquadr2 = (0d0,0d0)
        cquadr3 = (0d0,0d0)
        do i = 1,itheta
          ctmp1 = (0d0,0d0)
          ctmp2 = (0d0,0d0)
          ctmp3 = (0d0,0d0)
          do j = 1, iphi
            dist = (x0-xigrid(j,i))**2 + (y0-yigrid(j,i))**2
     $             + (z0-zigrid(j,i))**2
            dist = dsqrt(dist)
            ctmp = 1/(dist**3) * sphfgrid(j,i) * wsigrid(j,i)
            ctmp1 = ctmp1 + (x0-xigrid(j,i)) * ctmp
            ctmp2 = ctmp2 + (y0-yigrid(j,i)) * ctmp
            ctmp3 = ctmp3 + (z0-zigrid(j,i)) * ctmp
  !         write(931,*) j,i, "dist,tmp_inner,ctmp"
  !         write(931,*) dist,tmp_inner,ctmp
          enddo
          cquadr1 = ctmp1 * whts_i(i) + cquadr1
          cquadr2 = ctmp2 * whts_i(i) + cquadr2
          cquadr3 = ctmp3 * whts_i(i) + cquadr3
        enddo    
        cquadr1 = cquadr1/(-2*iphi)
        cquadr2 = cquadr2/(-2*iphi)
        cquadr3 = cquadr3/(-2*iphi)
              
      end SUBROUTINE       
     
