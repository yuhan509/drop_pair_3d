      subroutine integral_eval_ufield(itheta,iphi,
     $ xigrid, yigrid, zigrid, wsigrid, fxigrid, fyigrid, fzigrid, 
     $ whtsi, x0, y0, z0, quadr1, quadr2, quadr3)

      implicit real *8 (a-h,o-z)
      integer :: itheta, iphi
      real *8 :: x0, y0, z0, res1, res2, res3
      real *8 :: xigrid(iphi,itheta), yigrid(iphi,itheta),
     $ zigrid(iphi,itheta), wsigrid(iphi,itheta),
     $ fxigrid(iphi,itheta), fyigrid(iphi,itheta), fzigrid(iphi,itheta)
      real *8 :: whtsi(itheta)

              quadr1 = 0d0
              quadr2 = 0d0
              quadr3 = 0d0        
              do i = 1,itheta
                tmp1 = 0d0
                tmp2 = 0d0
                tmp3 = 0d0
                do j = 1, iphi
                  dist = (x0-xigrid(j,i))**2 
     $             + (y0-yigrid(j,i))**2
     $             + (z0-zigrid(j,i))**2
                  dist = dsqrt(dist)
              tmp = (x0-xigrid(j,i)) * fxigrid(j,i)
     $ + (y0-yigrid(j,i))* fyigrid(j,i)
     $ + (z0-zigrid(j,i))* fzigrid(j,i)
              tmp1 = tmp1 + (1d0/dist * fxigrid(j,i)
     $ + (x0-xigrid(j,i))/dist**3 * tmp) *wsigrid(j,i)
              tmp2 = tmp2 + (1d0/dist * fyigrid(j,i)
     $ + (y0-yigrid(j,i))/dist**3 * tmp)*wsigrid(j,i) 
              tmp3 = tmp3 + (1d0/dist * fzigrid(j,i)
     $ + (z0-zigrid(j,i))/dist**3 * tmp) *wsigrid(j,i)  
!            write(941,*) i,j, "dist, tmp, tmp1, tmp2, tmp3"
!           write(941,*)  dist, tmp, tmp1, tmp2, tmp3
                enddo
                quadr1 = tmp1 * whtsi(i) + quadr1
                quadr2 = tmp2 * whtsi(i) + quadr2
                quadr3 = tmp3 * whtsi(i) + quadr3
              enddo    
              quadr1 = quadr1/(4*iphi)
              quadr2 = quadr2/(4*iphi)
              quadr3 = quadr3/(4*iphi)
      end subroutine
c
c
c
c
      SUBROUTINE nearSing_integral_ufield(
     $ iterms, itheta, iphi,
     $ xigrid, yigrid, zigrid, wsigrid, 
     $ fxigrid, fyigrid, fzigrid, whtsi,
     $ umpole, dist_cl, sphf,  !! interpolate singluar integral
     $ xs, ys, zs, x0, y0, z0, quadr1, quadr2, quadr3)
         
         use mod_orgHelper     !! need to change name: ctheta... 
   
         implicit real *8 (a-h,o-z)
         integer :: iterms, itheta, iphi
         real *8 :: x0, y0, z0, xs, ys, zs, quadr1, quadr2, quadr3
         !! needed in integrand
         real *8 :: xigrid(iphi,itheta), yigrid(iphi,itheta),
     $ zigrid(iphi,itheta), wsigrid(iphi,itheta),
     $ fxigrid(iphi,itheta), fyigrid(iphi,itheta), fzigrid(iphi,itheta)
    
         real *8 :: whtsi(itheta)
         !! needed in singular integral interpolation
         complex *16 :: umpole(0:nterms, 0:nterms, 3), 
     $ sphf(0:nterms, -nterms:nterms)
   
         real *8 :: dist_cl
         real *8, ALLOCATABLE :: xarray(:)
         complex * 16, ALLOCATABLE :: farray(:,:)
         allocate(xarray(nintpol),farray(nintpol,3))

      !   PI = 4.D0*DATAN(1.D0)
      !   h = PI/(nterms + 1)/20

         if (dist_cl .ge. h) then
           call integral_eval_ufield(itheta,iphi,
     $ xigrid, yigrid, zigrid, wsigrid, fxigrid, fyigrid, fzigrid, 
     $ whtsi, x0, y0, z0, quadr1, quadr2, quadr3)

         else
         
         !! interpolate singluar integral from expansion
         tmp1 = 0d0
         tmp2 = 0d0
         tmp3 = 0d0
         do n = 0, nterms
            m = 0
            tmp1 = tmp1 + dble(umpole(n, m, 1) * sphf(n, m))
            tmp2 = tmp2 + dble(umpole(n, m, 2) * sphf(n, m))
            tmp3 = tmp3 + dble(umpole(n, m, 3) * sphf(n, m))
           do m = 1, n
             tmp1 = tmp1 + 2*dble(umpole(n, m, 1) * sphf(n, m))
             tmp2 = tmp2 + 2*dble(umpole(n, m, 2) * sphf(n, m))
             tmp3 = tmp3 + 2*dble(umpole(n, m, 3) * sphf(n, m))
           enddo
         enddo
         xarray(1) = 0d0 !! point on the surface
         farray(1,1) = tmp1
         farray(1,2) = tmp2
         farray(1,3) = tmp3
         
         ratio_x = (x0 - xs) / dist_cl
         ratio_y = (y0 - ys) / dist_cl
         ratio_z = (z0 - zs) / dist_cl
   
         !! call _eval for other regular/upsamp integrals
         PI = 4.D0*DATAN(1.D0)
         h = PI/(nterms + 1)
 !        write(903,*) "h=", h
         do i = 2, nintpol 
           d = dsqrt(dble(i-1)) * hintpol   !/2 !! modified  
           xarray(i) = d + xarray(1)
           x = xs + ratio_x * d
           y = ys + ratio_y * d
           z = zs + ratio_z * d
           call integral_eval_ufield(itheta,iphi,
     $ xigrid, yigrid, zigrid, wsigrid, fxigrid, fyigrid, fzigrid, 
     $ whtsi, x, y, z, quadr1, quadr2, quadr3)
           farray(i,1) = quadr1
           farray(i,2) = quadr2
           farray(i,3) = quadr3
         enddo

         !! target point
         xi = dist_cl
   
      call lagrange_value_1d(nintpol,xarray,farray(:,1), 1, xi, quadr1)
      call lagrange_value_1d(nintpol,xarray,farray(:,2), 1, xi, quadr2)
      call lagrange_value_1d(nintpol,xarray,farray(:,3), 1, xi, quadr3) 

        endif

          END SUBROUTINE
