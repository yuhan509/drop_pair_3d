cccccc Evaluate spherical harmonic function's on shifted grid in the 
cccccc pseudo time marching in the reparametrization.
cccccc Not a fast method.
c
c  Our definition of complex spherical harmonics is
c
c  Ynm(theta,phi)= sqrt(2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(im phi), 
c  Yn,-m(theta,phi) = sqrt(2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(-im phi),   for m >= 0.
c       
c  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.



      SUBROUTINE sph_grid(nterms, nphi, ntheta, thetaGrid, phiGrid,
     $ sphGrid) 
       implicit real *8 (a-h,o-z)
       integer :: nterms, nphi, ntheta
       complex *16 :: sphGrid(nphi,ntheta,0:nterms,-nterms:nterms)
       real *8 :: thetaGrid(nphi, ntheta), phiGrid(nphi, ntheta) 
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

c$OMP PARALLEL 
c$OMP do schedule(static) private(i,n,m,phij,ynm)      
      do j = 1, nphi
       do i = 1, ntheta
          phij = phiGrid(j,i)
          !write(47,*) "ij",i,j, phiGrid(j,i)
          call ylgndrf(nterms,dcos(thetaGrid(j,i)),ynm,rat1,rat2)
          do n = 0, nterms
            do m = -n, n
              !write(47,*) "ijnm",i,j,n,m, phiGrid(j,i)
              sphGrid(j,i,n,m) = ynm(n,abs(m))
     $         *DCMPLX(dcos(m*phij),dsin(m*phij))
            enddo
          enddo
       enddo
      enddo       
c$OMP end do
c$OMP END PARALLEL        

      END SUBROUTINE
