
ccc   Build spherical harmonic function on the grid
ccc   and all rotated sph function grids   
      SUBROUTINE rota_sphf(nterms, ntheta, nphi, nbeta, beta, nrot,
     $ ctheta, whts, ynms, zwsave, rotmat, 
     $ sphfgrid, rsphfgrids) 
cc input:
cc nterms, ntheta, nphi, nbeta, beta, nrot, ctheta, ynms, zwsave, rotmat
cc output:
cc sphfgrid, rsphfgrids     
      implicit real *8 (a-h,o-z)
      integer :: nterms, ntheta, nphi, nbeta, nrot
      real *8 ::  beta(nbeta), ctheta(ntheta), whts(ntheta)
      real *8 :: ynms(0:nterms,0:nterms,ntheta/2+1)
      complex *16 :: zwsave(4*nphi+15),
     $ rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta)
      complex *16 sphfgrid(nphi,ntheta,0:nterms,0:nterms),
     $ rsphfgrids(nphi,ntheta,nrot,nbeta,0:nterms,0:nterms)
      
      complex *16, ALLOCATABLE :: mpole(:,:)

      allocate(mpole(0:nterms,-nterms:nterms))
      
      call rotgrid_fsr_cmpl_init(nterms,nbeta,beta,rotmat)
            
      do n = 0, nterms
        do m = -n, n      
          mpole(n,m) = (0d0, 0d0)
        enddo
      enddo

    !  write(*,*) "beta(2)",beta(2),"rotmat(2,-1,-1,2)",rotmat(2,-1,-1,2)
      do n = 0, nterms
        do m = 0, n
          mpole(n,m) = (1d0, 0d0)
          write(*,*) "start building sph", n, m
          call sphtrans_cmpl(nterms,mpole,nphi,ntheta,sphfgrid(1,1,n,m),
     $     ctheta,ynms,zwsave)
           write(*,*) "start rotating sph", n, m
          call rotgrid_fsr_cmpl(nterms, mpole, nphi, ntheta,
     $        nbeta, beta, nrot, rsphfgrids(1,1,1,1,n,m),
     $        rotmat,ctheta,ynms,zwsave)             
          mpole(n,m) = (0d0, 0d0)
        enddo
      enddo
      write(*,*) "rotate spherical harmonics done!"
    !  write(*,*) nphi,ntheta,'ynms',ynms(3,2,1), ctheta(3), zwsave(10)
    !  write(*,*) sphfgrid(1,2,3,3),rsphfgrids(1,2,2,2,2,2)      
      
      END SUBROUTINE
c
c
c      
ccc   Build spherical harmonic function on the grid
      SUBROUTINE grid_sphf(nterms, ntheta, nphi,
     $ ctheta, ynms, zwsave, sphfgrid) 
cc input:
cc nterms, ntheta, nphi, ctheta, ynms, zwsave
cc output:
cc sphfgrid
      implicit real *8 (a-h,o-z)
      integer :: nterms, ntheta, nphi
      real *8 :: ctheta(ntheta)
      real *8 :: ynms(0:nterms,0:nterms,ntheta/2+1)
      complex *16 :: zwsave(4*nphi+15)
      complex *16 :: sphfgrid(nphi,ntheta,0:nterms,-nterms:nterms)
      
      complex *16, ALLOCATABLE :: mpole(:,:)

      allocate(mpole(0:nterms,-nterms:nterms))
      
      do n = 0, nterms
        do m = -n, n
          mpole(n,m) = (1d0, 0d0)
         call sphtrans_cmpl(nterms,mpole,nphi,ntheta,sphfgrid(:,:,n,m),
     $     ctheta,ynms,zwsave)      
          mpole(n,m) = (0d0, 0d0)
        enddo
      enddo

      write(*,*) "build spherical harmonics done!"
    !  write(*,*) nphi,ntheta,'ynms',ynms(1,2,1), ctheta(2), zwsave(10)
    !  write(*,*) sphfgrid(1,2,2,2)
      END SUBROUTINE
