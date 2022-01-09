        subroutine rotgrid_fsr_real_mpolesout(nterms,mpole,nphi,ntheta,
     $     nbeta,beta,nrot,grids,rotmat,ctheta,ynms,wsave,rfmpoles)
        implicit real *8 (a-h,o-z)
c
c  Rotate the real spherical harmonic grids.
c
c  Fast, FFT-based algorithm for rotating real spherical harmonic grids
c  (dimensioned NPHI-by-NTHETA) into new pole locations (beta_j, alpha_k), 
c  where alpha_k = 2*pi * k/nrot, k=0..nrot-1,  beta_j, j=1..nbeta. 
c
c  GRIDS = rotgrid_fsr_real(NTERMS,MPOLE,NPHI,NTHETA,NBETA,BETA,NROT,...
c       ROTMAT,CTHETA,YNMS,WSAVE) rotates 
c  the real spherical harmonics expansion of degree NTERMS 
c  about the z-axis by degree ALPHA_K and about the y-axis by degree BETA.
c  into a collection of new pole locations (beta_j, alpha_k) in spherical 
c  coordinates (theta, phi), where alpha_k = 2*pi * k/nrot, k=0..nrot-1,
c  beta_j, j=1..nbeta.
c
c  The rotated poles form a uniformly spaced grid on lattitude \theta.
c
c  grids - function values on the rotated grids, 
c             NPHI-by-NTHETA-by-NROT-by-NBETA real*8 matrix
c
c      Input parameters:
c
c  nterms - the number of terms in spherical harmonics expansion
c  mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,0:nterms)
c  nphi - the number of points in latitude discretization (for spherical grid)
c  ntheta - the number of points in meridian discretization (for spherical grid)
c  ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c  ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c  wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c
c  nbeta - the number of points in pole location meridian discretization
c  beta - angles for new pole locations beta_j, j=1..nbeta
c  nrot - angles for new pole locations alpha_k = 2*pi * k/nrot, k=0..nrot-1.
c       
c  rotmat - The rotation operators for directions (beta_j,0). 
c 
c  [rotmat, ctheta, ynms, wsave] must be initialized via 
c              a preceding call to rotgrid_fsr_real_init
c
c      Output parameters:
c
c  grids - function values on the rotated grids, 
c             NPHI-by-NTHETA-by-NROT-by-NBETA real*8 matrix
c
c  Our definition of complex spherical harmonics is
c
c  Ynm(theta,phi)= sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(im phi), 
c  Yn,-m(theta,phi) = sqrt( 2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(-im phi),   for m >= 0.
c       
c  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.
c
        complex *16 mpole(0:nterms,0:nterms)
        complex *16 rfmpoles(0:nterms,0:nterms,nrot,nbeta)

        real *8 beta(nbeta)
        real *8 rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta) 
        real *8 ctheta(ntheta),ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
c
        real *8 grids(nphi,ntheta,nrot,nbeta)
c
        complex *16, allocatable :: marray(:,:,:) 
c
        allocate( marray(0:nterms,0:nterms,nrot) )
c
        do k=1,nbeta
        call rot1lat_wfft_real
     $     (beta(k),nrot,nterms,
     $     nterms,nterms,mpole,nterms,marray,
     $     nterms,rotmat(0,-nterms,-nterms,k),nterms)
        do i=1,nrot
        call sphtrans_real(
     $     nterms,marray(0,0,i),
     $     nphi,ntheta,grids(1,1,i,k),
     $     ctheta,ynms,wsave)
           do n = 0, nterms
            do m = 0, n
                rfmpoles(n,m,i,k) = marray(n,m,i)
            enddo
           enddo
        enddo
        enddo
c
        return
        end