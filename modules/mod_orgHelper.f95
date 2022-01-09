MODULE mod_orgHelper
 use mod_assign

 real* 8, ALLOCATABLE:: ctheta(:), whts(:), theta(:), ynms(:,:,:) 
 complex *16, ALLOCATABLE::  zwsave(:), dwsave(:)
 real* 8, ALLOCATABLE:: swhts(:)
 SAVE
 
 contains
 
 SUBROUTINE mod_orgHelper_init()
  implicit real *8 (a-h,o-z)
      allocate(ctheta(ntheta), whts(ntheta), theta(ntheta))
      allocate(ynms(0:nterms,0:nterms,ntheta/2+1))
      allocate(zwsave(4*nphi+15))
      allocate(dwsave(4*nphi+15))
      allocate(swhts(ntheta))
      
  call grid_lege_init(nterms, ntheta, nphi, ctheta, theta, &
  whts, ynms, dwsave, zwsave)     
  
!cccccccccccccccccccccccccccccccccccccccccccccc
!cc    special gauss-legendre quadrature weights
!cccccccccccccccccccccccccccccccccccccccccccccc      
      do i = 1, ntheta
        swhts(i) = 0d0
        do n = 0, nterms
          call legepol(ctheta(i),n,pol,der)
!          write(*,*) 'x=',ctheta(j),'n=',i,'Pn(x)=',pol
          swhts(i) = swhts(i) + pol !/dsqrt(4d0*i + 2d0)
        enddo          
        swhts(i) = swhts(i) * whts(i) * 2d0 * dsin(theta(i)/2d0)
      enddo    
 
 END SUBROUTINE
      
      
END MODULE 
