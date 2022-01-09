      subroutine repara_real(nqterms, nqtheta, nqphi, idrop, istep,
     $  nup, thetagrid, phigrid,
     $ xqmpole, yqmpole, zqmpole, nextGrid)

      use mod_upsamp
      use OMP_LIB

      implicit real *8 (a-h,o-z)    
      integer :: nqterms, nqtheta, nqphi, nup
      real *8 :: thetagrid(nqphi,nqtheta), phigrid(nqphi,nqtheta)

      complex *16 :: xqmpole(0:nqterms,0:nqterms), 
     $ yqmpole(0:nqterms,0:nqterms), zqmpole(0:nqterms,0:nqterms)
      real *8 :: nextGrid(nqphi,nqtheta,3)
    
      complex *16, ALLOCATABLE :: xtmpole(:,:),ytmpole(:,:),ztmpole(:,:) 
      real *8, ALLOCATABLE :: ctheta(:), whts(:), ynms(:,:,:)
      complex *16, ALLOCATABLE :: dwsave(:) 
      real* 8, ALLOCATABLE :: dthdtau(:,:), dphdtau(:,:) 


      allocate(ctheta(nqtheta), whts(nqtheta))
      allocate(ynms(0:nqterms,0:nqterms,nqtheta/2+1))
      allocate(dwsave(4*nqphi+15))
      allocate(xtmpole(0:nqterms,0:nqterms))
      allocate(ytmpole(0:nqterms,0:nqterms))
      allocate(ztmpole(0:nqterms,0:nqterms))
      allocate(dthdtau(nqphi,nqtheta), dphdtau(nqphi,nqtheta))
      
      START = omp_get_wtime()
      do n = 0, nqterms
        do m = 0, n
          xtmpole(n,m) = xqmpole(n,m)
          ytmpole(n,m) = yqmpole(n,m)
          ztmpole(n,m) = zqmpole(n,m)
        enddo
      enddo

      !! "loaded ctheta, whts, ynms, dwsave"
      do i = 1, nqtheta
        ctheta(i) = up_ctheta(i,nup) 
        whts(i) = up_whts(i,nup)
      enddo
      do j = 1, nqphi*4+15
        dwsave(j) = up_dwsave(j,nup)
      enddo
      do n = 0, nqterms
        do m = 0, n
          do i = 1, nqtheta/2+1
             ynms(n,m,i) = up_ynms(n,m,i,nup)
          enddo
        enddo
      enddo 


      !!write(*,*) "loaded ctheta, whts, ynms, dwsave"

      !! do loop for repara_march
      iwrite = idrop * 100 + 50
      write(iwrite, *) "istep = ", istep

      do k = 1, i_repara_max
       !!write(*,*) "k = ", k

      !! write(*,*) " call eval_n_cutoff_real"

      !!! test n_cutoff = 0 diable : eval_n_cutoff_real
      call eval_n_cutoff_real(xtmpole, ytmpole, ztmpole, nqterms,
     $ ratio_cutoff, n_cutoff)
       write(iwrite, *) k, "n_cutoff=", n_cutoff

       !!write(*,*) " call eval_para_vel"



cccc  call eval_para_vel
      !!! debug
      iwrite = idrop*1000 + 501
      write(iwrite,*) "istep=",istep, "k=",k
      !!! test n_cutoff = 0
      ! n_cutoff = 0
      !!!!
       call eval_para_vel_real(nqterms, nqtheta, nqphi,
     $ idrop, n_cutoff, nup,
     $ xtmpole, ytmpole, ztmpole,
     $ dthdtau, dphdtau)
!      iwrite = idrop * 100 + 70
!      write(iwrite, *) "istep =",istep,"i_repara",k,"dthdtau, dphdtau"
!      test_cut = 1d-12
!      do i = 1, nqtheta
!        do j = 1, nqphi 
!          write(iwrite, *) i,j, dthdtau(j,i), dphdtau(j,i)
!          tmp = mod(dthdtau(j,i),test_cut)
!          dthdtau(j,i) = dthdtau(j,i) - mod(dthdtau(j,i),test_cut)
!          dphdtau(j,i) = dphdtau(j,i) - mod(dphdtau(j,i),test_cut)
!          !! forced symmetry
!          dthdtau(j,i) = dthdtau(1,i)
!          dphdtau(j,i) = 0d0
!        enddo 
!      enddo 
!!      write(iwrite, *) "cutoff less than", test_cut
!!      do i = 1, nqtheta
!!        do j = 1, nqphi 
!!          write(iwrite, *) i,j, dthdtau(j,i), dphdtau(j,i)
!!        enddo 
!!      enddo 

      !!write(*,*) " call repara_march_real"
      !!! need to record the last thetagrid and phigrid


      call repara_march_real(dtau, nqterms, nqtheta, nqphi, 
     $ thetagrid, phigrid,
     $ dthdtau, dphdtau, xqmpole, yqmpole, zqmpole, nextGrid)

!!       iwrite = idrop * 100 + 73
!!      write(iwrite, *) "istep =",istep,"i_repara",k,"nextGrid x,y,z"
!!      iwrite = idrop * 100 + 74
!!      write(iwrite, *) "istep =",istep,"i_repara",k,"nextGrid z ave dev"
!!      test_cut = 1d-14
!!      pi = 4d0*atan(1d0)
!!      do i = 1, nqtheta
!!        iwrite = idrop * 100 + 73
!!        !! find average
!!        sumz = 0d0
!!        do j = 1, nqphi 
!!          write(iwrite, *) i,j, nextGrid(j,i,1), nextGrid(j,i,2)
!!     $ , nextGrid(j,i,3)!, pi/nqphi*(j-1) 
!!          sumz = sumz + nextGrid(j,i,3)
!!        enddo 
!!        avez = sumz / nqphi
!!        !! standard deviation
!!        devz = 0d0
!!        do j = 1, nqphi 
!!          devz = devz + (nextGrid(j,i,3) - avez) ** 2
!!        enddo  
!!        devz = dsqrt(devz / nqphi)  
!!
!!        iwrite = idrop * 100 + 74  
!!        write(iwrite, *) "itheta=",i, avez, devz   
!!       enddo 
     
cccc  update xt,yt,ztmpole by nextGrid and the original parameters
      !!write(*,*) "call update_xyzmpole_real"
      call update_xyzmpole_real(nqterms, nqtheta, nqphi, 
     $ ctheta, whts, ynms, dwsave, nextGrid,
     $ xtmpole, ytmpole, ztmpole)

!!      iwrite = idrop * 100 + 75  
!!      write(iwrite, *) "istep =",istep,"i_repara",k,"ztmpole"
!!      do n = 0, nqterms
!!        do m = 0, n
!!          write(iwrite,*) n,m,ztmpole(n,m)
!!        enddo
!!      enddo
      !!write(*,*) "call eval_grid_qualitymeasure_real"
cccc  call eval_qua_measure 
      call eval_grid_qualitymeasure_real(nqterms, n_cutoff,
     $ xtmpole,ytmpole,ztmpole,quaMes)
     
       iwrite = idrop * 100 + 50
       write(iwrite,*) "quaMes",quaMes 

      enddo 

      !! return the updated coeff
      do n = 0, nqterms
        do m = 0, n
          xqmpole(n,m) = xtmpole(n,m)
          yqmpole(n,m) = ytmpole(n,m)
          zqmpole(n,m) = ztmpole(n,m)
        enddo
      enddo

      END SUBROUTINE
c
c
c
      SUBROUTINE eval_n_cutoff_real(xmpole, ympole, zmpole, nqterms,
     $ ratio_cutoff, n_cutoff)
     
       implicit real *8 (a-h,o-z)
       complex *16 xmpole(0:nqterms,0:nqterms)
       complex *16 ympole(0:nqterms,0:nqterms)
       complex *16 zmpole(0:nqterms,0:nqterms)  
       real *8 ratio_cutoff
       integer nqterms, n_cutoff

       real *8, ALLOCATABLE :: xyzmod(:)
       allocate(xyzmod(0:nqterms))

       sum1 = 0d0
       do n = 1,nqterms !! start with n = 1   ??? try 2 ???

       xyzmod(n) = 0d0
       m = 0
       xyzmod(n)=xyzmod(n)+dsqrt(dble(xmpole(n,m)*dconjg(xmpole(n,m)) 
     $                  + ympole(n,m)*dconjg(ympole(n,m))
     $                  + zmpole(n,m)*dconjg(zmpole(n,m))))
        do m = 1, n
       xyzmod(n)=xyzmod(n)+2*dsqrt(dble(xmpole(n,m)*dconjg(xmpole(n,m)) 
     $                  + ympole(n,m)*dconjg(ympole(n,m))
     $                  + zmpole(n,m)*dconjg(zmpole(n,m))))
        enddo
        sum1 = sum1 + xyzmod(n)
       enddo

       sum2 = 0d0
       n = nqterms !! start with n = 1
       thresh = sum1 * ratio_cutoff
       do while(sum2 .le. thresh) 
         sum2 = sum2 + xyzmod(n)         
         n = n - 1
       ENDDO

       n_cutoff = n + 2   

       if (n_cutoff .gt. nqterms) then
         write(*,*) "ratio_cutoff= ", ratio_cutoff,
     $ "is too small, set n_cutoff = nqterms"
         n_cutoff = nqterms
       endif
       return       

      END SUBROUTINE
c
c
c
      SUBROUTINE eval_n_cutoff_real_test(nqterms,ratio,arr,n_cut)
      implicit real *8 (a-h,o-z)
      integer :: nqterms,n_cut
      real *8 arr(nqterms), ratio

      sum1 = 0d0
      do n = 1, nqterms
        sum1 = sum1 + arr(n)
      enddo

      n = nqterms
      sum2 = 0d0
      do while(sum2 .le. sum1 * ratio) 
        sum2 = sum2 + arr(n)         
        n = n - 1
      ENDDO
       n_cut = n + 2    

       if (n_cut .gt. nqterms) then
            write(*,*) "ratio_cutoff is too small, set n_cut = nqterms"
            n_cut = nqterms
       endif

       return     
      END SUBROUTINE
c
c
c
      SUBROUTINE eval_grid_qualitymeasure_real(nqterms,n_cutoff,
     $ xmpole,ympole,zmpole,quaMes)
      implicit real *8 (a-h,o-z)
      integer :: nqterms,npterms,n_cutoff
      complex *16 xmpole(0:nqterms,0:nqterms)
      complex *16 ympole(0:nqterms,0:nqterms)
      complex *16 zmpole(0:nqterms,0:nqterms)  
      real *8 :: quaMes
      
        quaMes = 0d0            
        do n = n_cutoff,nqterms  !! perfect low-pass filter
          m = 0
          quaMes = quaMes
     $     + dble(xmpole(n,m)*dconjg(xmpole(n,m)))             
     $     + dble(ympole(n,m)*dconjg(ympole(n,m))) 
     $     + dble(zmpole(n,m)*dconjg(zmpole(n,m)))          
          do m = 1, n
           quaMes = quaMes
     $      + 2 * dble(xmpole(n,m)*dconjg(xmpole(n,m)))             
     $      + 2 * dble(ympole(n,m)*dconjg(ympole(n,m))) 
     $      + 2 * dble(zmpole(n,m)*dconjg(zmpole(n,m)))         
          enddo
        enddo 
        !write(50,*) "quaMes",quaMes                
      END SUBROUTINE
c
c
c
      SUBROUTINE update_xyzmpole_real(nqterms, nqtheta, nqphi,
     $ ctheta, whts, ynms, dwsave, nextGrid,
     $ xmpole, ympole, zmpole)

      integer :: nqterms, nqtheta, nqphi
      real* 8 :: ctheta(nqtheta), whts(nqtheta)
      real* 8 :: ynms(0:nqterms,0:nqterms,nqtheta/2+1)
      complex *16 :: dwsave(4*nqphi+15)      
      complex *16 xmpole(0:nqterms,0:nqterms)
      complex *16 ympole(0:nqterms,0:nqterms)
      complex *16 zmpole(0:nqterms,0:nqterms) 
      real *8 :: nextGrid(nqphi,nqtheta,3)    
     
      call sphtrans_fwd_real_rsv(nqterms,xmpole,nqphi,nqtheta,
     $ nextGrid(:,:,1),ctheta,whts,ynms,dwsave)
      call sphtrans_fwd_real_rsv(nqterms,ympole,nqphi,nqtheta,
     $ nextGrid(:,:,2),ctheta,whts,ynms,dwsave)    
      call sphtrans_fwd_real_rsv(nqterms,zmpole,nqphi,nqtheta,
     $ nextGrid(:,:,3),ctheta,whts,ynms,dwsave) 
      
      END SUBROUTINE
