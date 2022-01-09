      SUBROUTINE para_march(dtau, nterms, ntheta, nphi, theta,
     $ xmpole, ympole, zmpole, xgrid, ygrid, zgrid) 
      implicit real *8 (a-h,o-z)
     
      integer :: nterms, ntheta, nphi
      real *8 :: theta(ntheta), dtau 
      !real *8 :: thetagrid(nphi, ntheta),phigrid(nphi, ntheta)
      complex *16 xmpole(0:nterms,-nterms:nterms)
      complex *16 ympole(0:nterms,-nterms:nterms)
      complex *16 zmpole(0:nterms,-nterms:nterms) 
      real *8, ALLOCATABLE :: nexttheta(:,:), nextphi(:,:)  
      real *8, ALLOCATABLE :: dthdtau(:,:),dphdtau(:,:)
      complex *16, ALLOCATABLE :: sphGrid(:,:,:,:)
      complex *16 :: xgrid(nphi, ntheta),
     $ ygrid(nphi, ntheta),zgrid(nphi, ntheta)
      
      allocate(dthdtau(nphi, ntheta),dphdtau(nphi, ntheta))
      allocate(nexttheta(nphi, ntheta), nextphi(nphi, ntheta))   
      allocate(sphGrid(nphi,ntheta,0:nterms,-nterms:nterms))
      
      PI = 4.D0*DATAN(1.D0)
      phi = 2*PI/nphi      
      do i = 1, ntheta
        do j = 1,nphi
        nexttheta(j,i) = theta(i)
        nextphi(j,i) = phi*(j-1)
        enddo
      enddo

      !!forward Euler Metheod
      rate = 1d-2
      tot_time = 3d0
      nstep = 10
      dtau = tot_time/nstep 
      do k = 1, nstep
      do i = 1, ntheta
       do j = 1,nphi       
       dthdtau(j,i) = rate*2d0*dcos(3d0*nextphi(j,i))
       dphdtau(j,i) = rate*3d0*dcos(3d0*nexttheta(j,i))
       nexttheta(j,i) = nexttheta(j,i) + dthdtau(j,i) * dtau
       nextphi(j,i) = nextphi(j,i) + dphdtau(j,i) * dtau
       !write(47,*) "ij",i,j, nexttheta(j,i),nextphi(j,i)
       enddo
      enddo
      enddo
      
      call sph_grid(nterms, nphi, ntheta, nexttheta, nextphi,
     $ sphGrid)
      
      do i = 1, ntheta
        do j = 1, nphi
          tmpx = 0d0
          tmpy = 0d0
          tmpz = 0d0
          do n = 0, nterms
            do m = -n, n
               tmpx = tmpx + dble(sphGrid(j,i,n,m)*xmpole(n,m))
               tmpy = tmpy + dble(sphGrid(j,i,n,m)*ympole(n,m))
               tmpz = tmpz + dble(sphGrid(j,i,n,m)*zmpole(n,m))
            enddo
          enddo
          xgrid(j,i) = dcmplx(tmpx,0d0)
          ygrid(j,i) = dcmplx(tmpy,0d0)
          zgrid(j,i) = dcmplx(tmpz,0d0) 
        enddo
      enddo
         
      END SUBROUTINE


      SUBROUTINE para_march_real(nterms, ntheta, nphi, theta,
     $ xmpole, ympole, zmpole, xgrid, ygrid, zgrid) 
      implicit real *8 (a-h,o-z)
     
      integer :: nterms, ntheta, nphi
      real *8 :: theta(ntheta), dtau 
      !real *8 :: thetagrid(nphi, ntheta),phigrid(nphi, ntheta)
      complex *16 xmpole(0:nterms,0:nterms)
      complex *16 ympole(0:nterms,0:nterms)
      complex *16 zmpole(0:nterms,0:nterms) 
      real *8, ALLOCATABLE :: nexttheta(:,:), nextphi(:,:)  
      real *8, ALLOCATABLE :: dthdtau(:,:),dphdtau(:,:)
      complex *16, ALLOCATABLE :: sphGrid(:,:,:,:)
      real *8 :: xgrid(nphi, ntheta),
     $ ygrid(nphi, ntheta),zgrid(nphi, ntheta)
      
      allocate(dthdtau(nphi, ntheta),dphdtau(nphi, ntheta))
      allocate(nexttheta(nphi, ntheta), nextphi(nphi, ntheta))   
      allocate(sphGrid(nphi,ntheta,0:nterms,-nterms:nterms))
      
      PI = 4.D0*DATAN(1.D0)
      phi = 2*PI/nphi      
      do i = 1, ntheta
        do j = 1,nphi
        nexttheta(j,i) = theta(i)
        nextphi(j,i) = phi*(j-1)
        enddo
      enddo

      !!forward Euler Metheod
      rate = 1d-2
      tot_time = 3d0
      nstep = 10
      dtau = tot_time/nstep 
      do k = 1, nstep
      do i = 1, ntheta
       do j = 1,nphi       
       dthdtau(j,i) = rate*2d0*dcos(3d0*nextphi(j,i))
       dphdtau(j,i) = rate*3d0*dcos(3d0*nexttheta(j,i))
       nexttheta(j,i) = nexttheta(j,i) + dthdtau(j,i) * dtau
       nextphi(j,i) = nextphi(j,i) + dphdtau(j,i) * dtau
       !write(47,*) "ij",i,j, nexttheta(j,i),nextphi(j,i)
       enddo
      enddo
      enddo
      
      call sph_grid(nterms, nphi, ntheta, nexttheta, nextphi,
     $ sphGrid)
      
      do i = 1, ntheta
        do j = 1, nphi
          tmpx = 0d0
          tmpy = 0d0
          tmpz = 0d0
          do n = 0, nterms
            m = 0
            tmpx = tmpx + dble(sphGrid(j,i,n,m) * xmpole(n,m))
            tmpy = tmpy + dble(sphGrid(j,i,n,m) * ympole(n,m))
            tmpz = tmpz + dble(sphGrid(j,i,n,m) * zmpole(n,m))           
            do m = 1, n
               tmpx = tmpx + 2 * dble(sphGrid(j,i,n,m) * xmpole(n,m))
               tmpy = tmpy + 2 * dble(sphGrid(j,i,n,m) * ympole(n,m))
               tmpz = tmpz + 2 * dble(sphGrid(j,i,n,m) * zmpole(n,m))
            enddo
          enddo
          xgrid(j,i) = tmpx
          ygrid(j,i) = tmpy
          zgrid(j,i) = tmpz 
        enddo
      enddo
         
      END SUBROUTINE
