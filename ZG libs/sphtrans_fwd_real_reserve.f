       subroutine sphtrans_fwd_real_rsv(nterms,mpole,nphi,ntheta,ggrid,
     $     ctheta,whts,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Real valued O(p^3) forward spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c       No aliasing (nphi must be >= 2*nterms+1)
c       This code will break if nphi < 2*nterms+1
c
c       Warning, fgrid is a copy of ggrid and destroyed by this routine!
c       Hence, ggrid is reserved 
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       ggrid - function values on the grid, NPHI-by-NTHETA real*8 matrix
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c
c
        complex *16 mpole(0:nterms,0:nterms)
        real *8 ggrid(nphi,ntheta)
c
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd,ima
        data ima/(0.0d0,1.0d0)/
c
        real *8, ALLOCATABLE:: fgrid(:,:)
c
        allocate(fgrid(nphi,ntheta))
cc      make a copy of ggrid in fgrid
        do i = 1,ntheta
         do j = 1,nphi
           fgrid(j,i) = ggrid(j,i)
         enddo
        enddo
c
c
        do n=0,nterms
        do m=0,n
        mpole(n,m)=0
        enddo
        enddo
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call dfftf(nphi,fgrid(1,i),wsave)
        enddo
c
c

        do k=1,ntheta/2
c
c       ... form FFTPACK compatible spherical grid
c
        kf=k
        kr=ntheta+1-k
c
        do n=0,nterms,2
        mpole(n,0)=mpole(n,0)+
     $     whts(k)*(fgrid(1,kf)+fgrid(1,kr))*ynms(n,0,k)
        do m=1,n,2
          mf=2*m
          cd=dcmplx(fgrid(mf,kf),fgrid(mf+1,kf))-
     $       dcmplx(fgrid(mf,kr),fgrid(mf+1,kr))
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
        enddo
        do m=2,n,2
          mf=2*m
          cd=dcmplx(fgrid(mf,kf),fgrid(mf+1,kf))+
     $       dcmplx(fgrid(mf,kr),fgrid(mf+1,kr))
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
        enddo
        enddo
c
        do n=1,nterms,2
        mpole(n,0)=mpole(n,0)+
     $     whts(k)*(fgrid(1,kf)-fgrid(1,kr))*ynms(n,0,k)
        do m=1,n,2
          mf=2*m
          cd=dcmplx(fgrid(mf,kf),fgrid(mf+1,kf))+
     $       dcmplx(fgrid(mf,kr),fgrid(mf+1,kr))
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
        enddo
        do m=2,n,2
          mf=2*m
          cd=dcmplx(fgrid(mf,kf),fgrid(mf+1,kf))-
     $       dcmplx(fgrid(mf,kr),fgrid(mf+1,kr))
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
        enddo
        enddo
c
        enddo
c
c
        if( mod(ntheta,2) .eq. 1 ) then
        k=ntheta/2+1
c
c       ... form FFTPACK compatible spherical grid
c
        do n=0,nterms
        mpole(n,0)=mpole(n,0)+whts(k)*fgrid(1,k)*ynms(n,0,k)
        do m=1,n
          mf=2*m
          mpole(n,+m)=mpole(n,+m)+
     $       whts(k)*dcmplx(fgrid(mf,k),fgrid(mf+1,k))*ynms(n,m,k)
        enddo
        enddo
c
        endif
c

ccc        return
c
c       ... multiply by quadrature weights
c
        scale=0.5d0/dble(nphi)
c
        do n=0,nterms
        do m=0,n
        mpole(n,m)=mpole(n,m)*scale
        enddo
        enddo
c
        deallocate(fgrid)
        return
        end



        subroutine sphtrans_fwd_cmpl_rsv(nterms,mpole,nphi,ntheta,ggrid,
     $     ctheta,whts,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) forward spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c       No aliasing (nphi should be >= 2*nterms+1)
c       This code will break if nphi < nterms+1
c
c       Warning, fgrid is destroyed by this routine!
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 ggrid(nphi,ntheta)
c
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd
c
c
        complex *16, ALLOCATABLE:: fgrid(:,:)
c
        allocate(fgrid(nphi,ntheta))
cc      make a copy of ggrid in fgrid
        do i = 1,ntheta
         do j = 1,nphi
           fgrid(j,i) = ggrid(j,i)
         enddo
        enddo

        do n=0,nterms
          do m=-n,n
            mpole(n,m)=0
          enddo
        enddo
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
            call zfftf(nphi,fgrid(1,i),wsave)
        enddo
c
c
        do k=1,ntheta/2

c
c       ... form FFTPACK compatible spherical grid
c
          kf=k
          kr=ntheta+1-k
c
          do n=0,nterms,2
            cd=fgrid(1,kf)+fgrid(1,kr)
            mpole(n,0)=mpole(n,0)+whts(k)*ynms(n,0,k)*cd
            do m=1,n,2
              mf=m+1
              cd=fgrid(mf,kf)-fgrid(mf,kr)
              mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
              mf=nphi-m+1
              cd=fgrid(mf,kf)-fgrid(mf,kr)
              mpole(n,-m)=mpole(n,-m)+whts(k)*ynms(n,m,k)*cd
            enddo
            do m=2,n,2
              mf=m+1
              cd=fgrid(mf,kf)+fgrid(mf,kr)
              mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
              mf=nphi-m+1
              cd=fgrid(mf,kf)+fgrid(mf,kr)
              mpole(n,-m)=mpole(n,-m)+whts(k)*ynms(n,m,k)*cd
            enddo
          enddo

          do n=1,nterms,2
            cd=fgrid(1,kf)-fgrid(1,kr)
            mpole(n,0)=mpole(n,0)+whts(k)*ynms(n,0,k)*cd
            do m=1,n,2
              mf=m+1
              cd=fgrid(mf,kf)+fgrid(mf,kr)
              mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
              mf=nphi-m+1
              cd=fgrid(mf,kf)+fgrid(mf,kr)
              mpole(n,-m)=mpole(n,-m)+whts(k)*ynms(n,m,k)*cd
            enddo
            do m=2,n,2
              mf=m+1
              cd=fgrid(mf,kf)-fgrid(mf,kr)
              mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
              mf=nphi-m+1
              cd=fgrid(mf,kf)-fgrid(mf,kr)
              mpole(n,-m)=mpole(n,-m)+whts(k)*ynms(n,m,k)*cd
            enddo
          enddo

        enddo
c
c
        if( mod(ntheta,2) .eq. 1 ) then
            k=ntheta/2+1
c
c       ... form FFTPACK compatible spherical grid
c
            do n=0,nterms
              mpole(n,0)=mpole(n,0)+whts(k)*fgrid(1,k)*ynms(n,0,k)
              do m=1,n
                mf=m+1
                mpole(n,+m)=mpole(n,+m)+whts(k)*fgrid(mf,k)*ynms(n,m,k)
                mf=nphi-m+1
                mpole(n,-m)=mpole(n,-m)+whts(k)*fgrid(mf,k)*ynms(n,m,k)
              enddo
            enddo
c
        endif
c
ccc        return
c
c       ... multiply by quadrature weights
c
        scale=0.5d0/dble(nphi)
c
        do n=0,nterms
            do m=-n,n
                mpole(n,m)=mpole(n,m)*scale
            enddo
        enddo
c
        deallocate(fgrid)
        return
        end
c       