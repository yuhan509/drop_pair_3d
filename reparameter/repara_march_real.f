      SUBROUTINE repara_march_real(dtau, nterms, ntheta, nphi, 
     $ thetagrid, phigrid, dthdtau, dphdtau, 
     $ xmpole, ympole, zmpole, nextGrid) 
      implicit real *8 (a-h,o-z)
     
      integer :: nterms, ntheta, nphi
      real *8 :: theta(ntheta), dtau 
      real *8 :: dthdtau(nphi, ntheta),dphdtau(nphi, ntheta)
      complex *16 xmpole(0:nterms,0:nterms)
      complex *16 ympole(0:nterms,0:nterms)
      complex *16 zmpole(0:nterms,0:nterms) 
      real *8 :: nextGrid(nphi, ntheta, 3)
      real *8:: thetagrid(nphi,ntheta), phigrid(nphi,ntheta)
      complex *16, ALLOCATABLE :: sphGrid(:,:,:,:)
      !complex *16:: sphharf(nphi,ntheta,0:nterms,-nterms:nterms)
      
      allocate(sphGrid(nphi,ntheta,0:nterms,-nterms:nterms))
      !!forward Euler Metheod
      do i = 1, ntheta
       do j = 1,nphi
       thetagrid(j,i) = thetagrid(j,i) + dthdtau(j,i) * dtau
       phigrid(j,i) = phigrid(j,i) + dphdtau(j,i) * dtau
       enddo
      enddo
      
      call sph_grid(nterms, nphi, ntheta, thetagrid, phigrid,
     $ sphGrid)

c$OMP PARALLEL 
c$OMP do schedule(static) private(i,n,m,tmpx,tmpy,tmpz)           
      do j = 1, nphi
        do i = 1, ntheta
          tmpx = 0d0
          tmpy = 0d0
          tmpz = 0d0
          do n = 0, nterms
            m = 0
            tmpx = tmpx + dble(sphGrid(j,i,n,m)*xmpole(n,m))
            tmpy = tmpy + dble(sphGrid(j,i,n,m)*ympole(n,m))
            tmpz = tmpz + dble(sphGrid(j,i,n,m)*zmpole(n,m))
            
            do m = 1, n
               tmpx = tmpx + 2 *dble(sphGrid(j,i,n,m)*xmpole(n,m))
               tmpy = tmpy + 2 *dble(sphGrid(j,i,n,m)*ympole(n,m))
               tmpz = tmpz + 2 *dble(sphGrid(j,i,n,m)*zmpole(n,m))
               !write(49,*) sphGrid(j,i,n,m) - sphharf(j,i,n,m)
            enddo
          enddo
          nextGrid(j,i,1) = tmpx
          nextGrid(j,i,2) = tmpy
          nextGrid(j,i,3) = tmpz 
        enddo
      enddo
c$OMP end do
c$OMP END PARALLEL                 
      END SUBROUTINE

