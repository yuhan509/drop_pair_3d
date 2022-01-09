      SUBROUTINE repara_real_upsamp( idrop, istep,
     $ xgrid, ygrid, zgrid, xmpole, ympole, zmpole)
     
      use mod_assign
      use mod_orgHelper
      use mod_upsamp
      use OMP_LIB

      implicit real *8 (a-h,o-z)    

      integer :: idrop, istep

      real *8:: xgrid(nphi,ntheta),ygrid(nphi,ntheta),zgrid(nphi,ntheta)
      complex *16 :: xmpole(0:nterms,0:nterms),
     $ ympole(0:nterms,0:nterms),zmpole(0:nterms,0:nterms)

      integer :: nqterms
      complex *16,ALLOCATABLE:: xqmpole(:,:), yqmpole(:,:), zqmpole(:,:)
      real* 8, ALLOCATABLE :: thetagrid(:,:), phigrid(:,:),
     $ nextGrid(:,:,:)   
      real* 8  :: START,END

      nqterms = nterms * nup_repara
      allocate(xqmpole(0:nqterms,0:nqterms))
      allocate(yqmpole(0:nqterms,0:nqterms))
      allocate(zqmpole(0:nqterms,0:nqterms))

      call grid_dim_init(nqterms, nqtheta, nqphi)
      allocate(thetagrid(nqphi,nqtheta), phigrid(nqphi,nqtheta))
      allocate(nextGrid(nqphi,nqtheta,3)) 

      !! x,y,zgrid fwd sph transform to x,y,zmpole
      call sphtrans_fwd_real_rsv(nterms,xmpole,nphi,ntheta,xgrid,
     $     ctheta,whts,ynms,dwsave)
      call sphtrans_fwd_real_rsv(nterms,ympole,nphi,ntheta,ygrid,
     $     ctheta,whts,ynms,dwsave)    
      call sphtrans_fwd_real_rsv(nterms,zmpole,nphi,ntheta,zgrid,
     $     ctheta,whts,ynms,dwsave)    
      
      !! test if x,y,zmpole can transform back to x,y,zgrid
      
      

      !! x,y,zmpole is the starting point of repara
      !! upsample x,y,zmpole to xq,yq,zqmpole
      iwrite = idrop*100 + 61
      write(iwrite, *) "istep = ", istep, "after marching before repara"
      do n = 0, nterms
       do m = 0, n
         xqmpole(n,m) = xmpole(n,m)
         yqmpole(n,m) = ympole(n,m)
         zqmpole(n,m) = zmpole(n,m)
         write(iwrite,*) n,m, xmpole(n,m), ympole(n,m), zmpole(n,m)
       enddo
      enddo
      
      do n = nterms + 1, nqterms
        do m = 0, n
            xqmpole(n,m) = (0d0,0d0)
            yqmpole(n,m) = (0d0,0d0)
            zqmpole(n,m) = (0d0,0d0)
        enddo
      enddo         
      !!write(*,*)  "upsamp to xq,yq,zqmpole"

      !! initialize the upsampled standard para grid thetagrid & phigrid
      PI = 4.D0*DATAN(1.D0)
      phi = 2*PI/nqphi      
      do i = 1, nqtheta
        do j = 1,nqphi
          thetagrid(j,i) = up_theta(i,nup_repara)
          phigrid(j,i) = phi*(j-1)
        enddo
      enddo

      !!write(*,*)  "start repara xq,yq,zqmpole"
      !! input :...
      !! output : xq,yq,zqmpole


      call repara_real(nqterms, nqtheta, nqphi, idrop, istep, 
     $ nup_repara, thetagrid, phigrid,
     $ xqmpole, yqmpole, zqmpole, nextGrid)

      !!write(*,*)  "repara to xq,yq,zqmpole"

      !! truncate the final xq, yq,zqmpole to x,y,zmpole
      iwrite = idrop*100 + 62
      write(iwrite, *) "istep = ", istep, "after repara"
      do n = 0, nterms
        do m = 0, n
          xmpole(n,m) = xqmpole(n,m)
          ympole(n,m) = yqmpole(n,m)
          zmpole(n,m) = zqmpole(n,m)
          write(iwrite,*) n,m, xmpole(n,m), ympole(n,m), zmpole(n,m)
          !! screening small values from x,y,zmpole after repara 
!!          tol_squ = 1d-30
!!          if (dble(xmpole(n,m))**2 .lt. tol_squ) then
!!            xmpole(n,m) = dcmplx(0d0, dimag(xmpole(n,m)))
!!          endif
!!          if (dble(ympole(n,m))**2 .lt. tol_squ) then
!!            ympole(n,m) = dcmplx(0d0, dimag(ympole(n,m)))
!!          endif
!!          if (dble(zmpole(n,m))**2 .lt. tol_squ) then
!!            zmpole(n,m) = dcmplx(0d0, dimag(zmpole(n,m)))
!!          endif
!!          if (dimag(xmpole(n,m))**2 .lt. tol_squ) then
!!            xmpole(n,m) = dcmplx(dble(xmpole(n,m)), 0d0)
!!          endif
!!          if (dimag(ympole(n,m))**2 .lt. tol_squ) then
!!            ympole(n,m) = dcmplx(dble(ympole(n,m)), 0d0)
!!          endif
!!          if (dimag(zmpole(n,m))**2 .lt. tol_squ) then
!!            zmpole(n,m) = dcmplx(dble(zmpole(n,m)), 0d0)
!!          endif
        enddo
      enddo 

!      do n = 0, nterms
!        do m = 0, n
!          write(iwrite,*) n,m, xmpole(n,m), ympole(n,m), zmpole(n,m)
!        enddo
!      enddo 


      !!
      END SUBROUTINE
