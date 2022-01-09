      SUBROUTINE eval_velo_single_wmod(istep, Ca, 
     $ x1grid, y1grid, z1grid, x1mpole, y1mpole, z1mpole, 
     $ h1grid, w1grid, n1grid, u1grids, en1mpole)
    
      use mod_orgHelper
      use mod_rotsphf
      use mod_upsamp
      
      implicit real *8 (a-h,o-z)
      
      integer :: istat = 2, istep
      real *8 :: Ca

      real *8 :: x1grid(nphi,ntheta),y1grid(nphi,ntheta),
     $ z1grid(nphi,ntheta),h1grid(nphi,ntheta),w1grid(nphi,ntheta), 
     $ n1grid(nphi,ntheta,3)
      complex *16 :: x1mpole(0:nterms,0:nterms),
     $ y1mpole(0:nterms,0:nterms),z1mpole(0:nterms,0:nterms)

      !! output 
      real *8 :: u1grids(nphi,ntheta,3)
      complex *16 :: en1mpole(0:nterms,-nterms:nterms)
      !! end of output
      
      complex *16, ALLOCATABLE :: sphf_s(:,:)
      integer :: iterms, itheta, iphi
      real *8, ALLOCATABLE :: whts_i(:)

      complex *16,dimension(:,:),ALLOCATABLE:: cEinfn1grid
      real *8, dimension(:,:), ALLOCATABLE:: ws1grid
      complex *16, dimension(:,:), ALLOCATABLE:: ws1mpole
      real *8, dimension(:,:,:,:), ALLOCATABLE:: rx1grids, ry1grids,
     $ rz1grids, rws1grids

      complex *16, dimension(:,:,:,:), ALLOCATABLE :: cls1d1grids
      complex *16, dimension(:,:,:,:), ALLOCATABLE :: cls1d1mpoles

      complex *16, dimension(:,:), ALLOCATABLE:: Einfn1mpole

      complex *16, ALLOCATABLE:: cefield_mat(:,:), cEinf_col(:) 

      INTEGER :: NRHS=1,nterms2
      INTEGER *4 :: INFO
      INTEGER, allocatable :: ipiv(:)

      real *8, ALLOCATABLE, dimension(:,:):: En1
      real *8, dimension(:,:), allocatable:: fx1,fy1,fz1
      complex *16, dimension(:,:), allocatable :: fx1mpole, fy1mpole,
     $ fz1mpole
      real *8, dimension(:,:,:,:), allocatable :: rfx1grids, rfy1grids,
     $ rfz1grids
      real *8, dimension(:,:,:), allocatable :: us1d1grid
      
      complex *16::ctmp,cquadr,ctmp1,ctmp2,ctmp3,cquadr1,cquadr2,cquadr3
     $ , cres

      !! testing
      complex *16, allocatable :: mpole(:,:),ggrid(:,:)
      !! dimension of matrix block
      nterms2=(nterms+1)*(nterms+1)
      
      PI = 4.D0*DATAN(1.D0)
cccccccccccccccccccccccccccccccccccccccccccccc
      allocate(ws1grid(nphi,ntheta))
      allocate(cEinfn1grid(nphi,ntheta))
      do i=1,ntheta
        do j=1,nphi
          ws1grid(j,i) = w1grid(j,i)/dsin(theta(i))
          !! E_inf * ngrid
          cEinfn1grid(j,i) = dcmplx(n1grid(j,i,3), 0d0)
          !write(88,*) j, i, ngrid(j,i,3), ctheta(i)
        enddo
      enddo

      allocate(ws1mpole(0:nterms,0:nterms))

      call sphtrans_fwd_real_rsv(nterms,ws1mpole,nphi,ntheta,ws1grid,
     $  ctheta,whts,ynms,dwsave)
cccccccccccccccccccccccccccccccccccc
cccccc rotate x,y,z mpole and local area element ccccccccccc
      nbeta = ntheta
      nrot = nphi

      allocate(rx1grids(nphi,ntheta,nrot,nbeta),
     $         ry1grids(nphi,ntheta,nrot,nbeta),
     $         rz1grids(nphi,ntheta,nrot,nbeta),
     $         rws1grids(nphi,ntheta,nrot,nbeta))

      call rotgrid_dsr_real(nterms,x1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rx1grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,y1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,ry1grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,z1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rz1grids,
     $ rotmat,ctheta,ynms,dwsave)     
      call rotgrid_dsr_real(nterms,ws1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rws1grids,
     $ rotmat,ctheta,ynms,dwsave)      

cc    special gauss-legendre quadrature(lapalce) for all rotated grids
      allocate(cls1d1grids(nrot,nbeta,0:nterms,-nterms:nterms))
      do n = 0,nterms
        do m = -n,n
          !!write(*,*) "using rotating sph", n, m
          do kk = 1,nrot
            do k = 1,nbeta  

              write(931,*) n,m,kk,k

              call integral_eval_efield_cmpl(nterms,ntheta,nphi,
     $        rx1grids(:,:,kk,k),ry1grids(:,:,kk,k),rz1grids(:,:,kk,k),
     $        rws1grids(:,:,kk,k), rsphfgrids(:,:,kk,k,n,m), swhts,
     $        x1grid(kk,k), y1grid(kk,k), z1grid(kk,k), n1grid(kk,k,:),
     $        cres)  
              !! naming : c- is for complex; l- is for laplace kernel;
              !! s1- means source(kernel) point is on drop 1;
              !! d1- means the contribution from integrating over drop 1.
              !! (kk,k) represents kernel point, will serve as 
              cls1d1grids(kk,k,n,m) = cres

              write(933,*) kk,k,n,m,cls1d1grids(kk,k,n,m)
!              write(13,*) 'n=',n,'m=',m,'irot=',kk,'jbeta=',k,
!     $  'special g-l quadr laplace',clapl_quadrs(kk,k,n,m)

              
            enddo
          enddo
        enddo
      enddo     

ccccccccccccccc


ccccccccc Inner Product ccccccccccccc
      allocate( mpole(0:nterms,-nterms:nterms), ggrid(nphi,ntheta),
     $ cls1d1mpoles(0:nterms,-nterms:nterms,0:nterms,-nterms:nterms))

      do n = 0,nterms
        do m = -n,n
          call sphtrans_fwd_cmpl_rsv(nterms,cls1d1mpoles(:,:,n,m),
     $    nphi,ntheta,cls1d1grids(:,:,n,m),ctheta,whts,ynms,zwsave)
        enddo
      enddo  

      allocate(Einfn1mpole(0:nterms,-nterms:nterms))
      call sphtrans_fwd_cmpl_rsv(nterms,Einfn1mpole,
     $    nphi,ntheta,cEinfn1grid,ctheta,whts,ynms,zwsave)

cccccccccccccccccccccccccccccccccccc

      allocate(cefield_mat(nterms2,nterms2))
      j = 1
      do n = 0,nterms !! start 0
        do m = -n, n
          i = 1
          do n1 = 0,nterms !! start 0
            do m1 = -n1, n1
              cefield_mat(i, j)=cls1d1mpoles(n1,m1,n,m)
!              write(15,*) 'n=',n,'m=',m,'n1=',n1,'m1=',m1,'i=',i,'j=',j,
!     $  cap_lpl(i,j)
              i = i + 1    
            enddo
          enddo
          j = j + 1
        enddo
      enddo    

      do i = 1, nterms2
        !! note: factor 4*pi is because the spherical harmonic functions we used are not normalized
         cefield_mat(i,i) = cefield_mat(i,i)+dcmplx(0.5d0,0d0)  
!         write(15, *) 'i=',i,cap_lpl(i,i)
      enddo

      allocate(cEinf_col(nterms2))
      write(14,*) 'Y proj of E_inf*ngrid, ALL B terms in A * X = B'
      write(21,*) "only nonzero B terms, in A * X = B"
      i = 1
      do n = 0, nterms !! start 0
          do m = -n, n              
              cEinf_col(i) = Einfn1mpole(n, m)
             ! write(14,*) 'n',n,'m',m,'i=',i,cap_Einf(i)
              i = i + 1
          enddo
      enddo

      do i = 1, nterms2
        if (abs(cEinf_col(i)*dconjg(cEinf_col(i))).lt.1d-25) then
          cEinf_col(i) = 0d0
        else
          write(21,*)  'n',n,'m',m,'i=',i, cEinf_col(i)
        endif
      enddo
      
cccccccccccc
       write(19,*)  " all A terms, in A * X = B"   
       write(20,*)  " only nonzero A terms, in A * X = B"  
       do k = 1, nterms2
         do l = 1, nterms2
           !write(19,*) 'i=',k,'j=',l, cap_lpl(k,l)
           tmp = abs(cefield_mat(k,l)*dconjg(cefield_mat(k,l)))
           if (tmp.lt.1d-25) then 
            cefield_mat(k,l) = 0d0
           else
            write(20,*) 'i=',k,'j=',l, cefield_mat(k,l)
           endif
          enddo
       enddo

      !!!! special case for the first step sphere shape
      !! avoid all zero entry on the first line of linear system (A and B) 
       if (istep .eq. 1) then 
        cefield_mat(1,1) = (1d0, 0d0)
      endif          
cccccccccccc
      allocate(ipiv(nterms2))         !! n start with 0            
      CALL zgesv(nterms2,NRHS,cefield_mat,nterms2,IPIV,
     $ cEinf_col,nterms2,INFO)      
!      CALL zgesv(N2-1,NRHS,cap_lpl,N2-1,
!     $ IPIV,cap_Einf,N2-1,INFO)      
      write(16,*) 'linear solver zgesv solution:','INFO=',INFO

      !!!! special case for the first step sphere shape
      if (istep .eq. 1) then 
        en1mpole(0,0) = (0d0, 0d0)
      endif      

      i = 1
      do n = 0, nterms !! start 0
          do m = -n, n              
              en1mpole(n, m) = cEinf_col(i)   
              write(16,*) 'i=',i,'n',n,'m',m,cEinf_col(i)
              tmp = abs(en1mpole(n, m)*dconjg(en1mpole(n, m)))
              if (tmp.lt.1d-25) then
                en1mpole(n, m) = 0d0
              endif
              i = i + 1
          enddo
      enddo

!      !!!! special case for the first step sphere shape
!      if (istep .eq. 1) then 
!        enmpole(0,0) = (0d0, 0d0)
!      endif
      
ccccc spherical harmonic transform ccccccc
      allocate(En1(nphi, ntheta))

      do i = 1,ntheta
          do j = 1, nphi
            ctmp1 = (0d0,0d0)
            do n = 0, nterms   !! start 0, only take odd terms
              do m = -n, n              
                  ctmp1 = ctmp1 + en1mpole(n, m)*sphfgrid(j,i,n,m)              
              enddo
            enddo 
            En1(j,i) = dble(ctmp1)
            write(17,*) 'j=',j,'i=',i, En1(j,i),En1(j,i)/ctheta(i)
          enddo
      enddo

cccccccccc total force mpole cccccccccc
c          ( En^2/2*Ca + 2*H ) * ngrid
ccccccccccccccccccccccccccccccccccccccc    

      allocate(fx1(nphi,ntheta), fy1(nphi,ntheta), fz1(nphi,ntheta))

      Caflow = 0.0d0
      do i = 1, ntheta   
      do j = 1, nphi    

          tmp1 = (En1(j,i)**2/2d0*Ca + 2d0*h1grid(j,i))
          fx1(j,i) = tmp1 * n1grid(j,i,1)
          fy1(j,i) = tmp1 * n1grid(j,i,2)
          fz1(j,i) = tmp1 * n1grid(j,i,3)
        write(32,*) j,i,En1(j,i), h1grid(j,i), n1grid(j,i,1),
     $ n1grid(j,i,2), n1grid(j,i,3)
        enddo
      enddo

      write(947,*)  "j,i, fx(j,i), fy(j,i), fz(j,i)"
      do i = 1, ntheta
      do j = 1, nphi

          write(947,*)  j,i, fx1(j,i), fy1(j,i), fz1(j,i)
        enddo
      enddo

      allocate(fx1mpole(0:nterms,0:nterms), fy1mpole(0:nterms,0:nterms),
     $ fz1mpole(0:nterms,0:nterms))
      
      call sphtrans_fwd_real_rsv(nterms,fx1mpole,nphi,ntheta,fx1,
     $     ctheta,whts,ynms,dwsave)
      call sphtrans_fwd_real_rsv(nterms,fy1mpole,nphi,ntheta,fy1,
     $     ctheta,whts,ynms,dwsave)    
      call sphtrans_fwd_real_rsv(nterms,fz1mpole,nphi,ntheta,fz1,
     $     ctheta,whts,ynms,dwsave)           
cccccccccccccccccccccccccccccccccccccccc
      write(945,*)  "n,m, fxmpole(n,m), fympole(n,m), fzmpole(n,m)"
      do n = 0, nterms
        do m = 0,n
          write(945,*)  n,m, fx1mpole(n,m), fy1mpole(n,m), fz1mpole(n,m)
        enddo
      enddo
cccccc rotate force components mpole ccccccccccc

      allocate(rfx1grids(nphi,ntheta,nrot,nbeta),
     $ rfy1grids(nphi,ntheta,nrot,nbeta),
     $ rfz1grids(nphi,ntheta,nrot,nbeta))

      call rotgrid_dsr_real(nterms,fx1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rfx1grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,fy1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rfy1grids,
     $ rotmat,ctheta,ynms,dwsave)
      call rotgrid_dsr_real(nterms,fz1mpole,nphi,ntheta,
     $ nbeta,theta,nrot,rfz1grids,
     $ rotmat,ctheta,ynms,dwsave)  
     
      do j = 1, nphi
        do i = 1, ntheta
          do kk = 1, nrot
            do k = 1, nbeta
              write(943,*) j,i,kk,k
              write(943,*) rfx1grids(j,i,kk,k),rfy1grids(j,i,kk,k)
     $ ,rfz1grids(j,i,kk,k), rws1grids(j,i,kk,k)
            enddo
          enddo
        enddo
      enddo
cccccccccccccccccccccccccccccccccccccccccccccc
      
      allocate(us1d1grid(nphi,ntheta,3))

          do kk = 1,nphi
            do k = 1,ntheta
              write(941,*) kk, k
              call integral_eval_ufield(ntheta,nphi,
     $ rx1grids(:,:,kk,k), ry1grids(:,:,kk,k), rz1grids(:,:,kk,k), 
     $ rws1grids(:,:,kk,k), rfx1grids(:,:,kk,k), rfy1grids(:,:,kk,k), 
     $ rfz1grids(:,:,kk,k), swhts, 
     $ x1grid(kk,k), y1grid(kk,k), z1grid(kk,k), 
     $ quadr1, quadr2, quadr3)
              us1d1grid(kk,k,1) = quadr1 
              us1d1grid(kk,k,2) = quadr2
              us1d1grid(kk,k,3) = quadr3 
              write(34,*) kk, k, us1d1grid(kk,k,1), us1d1grid(kk,k,2),
     $  us1d1grid(kk,k,3)
            enddo
          enddo 
          
          do j = 1, nphi
            do i = 1, ntheta
              do k = 1,3
              u1grids(j,i,k) = us1d1grid(j,i,k) 
              enddo
            enddo
          enddo

      deallocate(cEinfn1grid)
      deallocate(ws1grid, ws1mpole,rx1grids, ry1grids,rz1grids, 
     $ rws1grids)
      deallocate(cls1d1grids)
      deallocate(cls1d1mpoles)
      deallocate(Einfn1mpole, cefield_mat, cEinf_col, ipiv)
      deallocate(En1, fx1, fy1, fz1)
      deallocate(fx1mpole, fy1mpole, fz1mpole)
      deallocate(rfx1grids, rfy1grids, rfz1grids)
      deallocate(us1d1grid)

      end SUBROUTINE

c
c  Our definition of complex spherical harmonics is
c
c  Ynm(theta,phi)= sqrt(2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(im phi), 
c  Yn,-m(theta,phi) = sqrt(2n+1) sqrt((n-m)!/(n+m)!) 
c                  Pnm(cos theta) e^(-im phi),   for m >= 0.
c       
c  Note that we do not include the Condon-Shortley phase (-1)^m, if m<0.




 
