module mod_eqsamp
 use mod_assign
 integer ::  nsamp 
 
 real *8, ALLOCATABLE,dimension(:,:) :: ethetagrid, ephigrid
 complex *16, ALLOCATABLE :: eyfgrid(:,:,:,:)
 real *8, ALLOCATABLE, dimension(:,:) :: exgrid, eygrid, ezgrid
 real *8, ALLOCATABLE, dimension(:,:) :: ex1grid, ey1grid, ez1grid
 real *8, ALLOCATABLE, dimension(:,:) :: ex2grid, ey2grid, ez2grid
 save

 contains
 
 SUBROUTINE mod_eqsamp_nsamp_init()
      implicit real *8 (a-h,o-z)
      if (mod(ntheta,2) .eq. 0) then
       nsamp = ntheta + 1
      else 
       nsamp = ntheta
      endif
      
      allocate(ephigrid(nphi,nsamp), ethetagrid(nphi,nsamp)) 
      allocate(eyfgrid(nphi,nsamp,0:nterms,-nterms:nterms))
      allocate(exgrid(nphi,nsamp), eygrid(nphi,nsamp), &
       ezgrid(nphi,nsamp))
      allocate(ex1grid(nphi,nsamp), ey1grid(nphi,nsamp), &
       ez1grid(nphi,nsamp))
      allocate(ex2grid(nphi,nsamp), ey2grid(nphi,nsamp), &
       ez2grid(nphi,nsamp))

      PI = 4.D0*DATAN(1.D0)
      phi = 2*PI/nphi 
      do i = 1, nsamp
       do j = 1, nphi
          ephigrid(j, i) = (j-1)*phi
          ethetagrid(j, i) = (i-1)*PI/(nsamp - 1)
       enddo
      enddo
      call sph_grid(nterms,nphi,nsamp,ethetagrid,ephigrid,eyfgrid)
      
 END SUBROUTINE 


 SUBROUTINE mod_eqsamp_eval_grid_xyz(istep, xmpole, ympole, zmpole)
  implicit real *8 (a-h,o-z)     
  integer :: istep
  real *8 :: longaxis, shortaxis
  complex *16 :: xmpole(0:nterms, 0:nterms),ympole(0:nterms, 0:nterms),&
  zmpole(0:nterms, 0:nterms)       
  
     do i = 1, nsamp
       do j = 1, nphi
        exgrid(j,i) = 0d0
        eygrid(j,i) = 0d0
        ezgrid(j,i) = 0d0
        do n = 0,nterms
         m = 0
         exgrid(j,i) = exgrid(j,i)   & 
           + dble(xmpole(n,m) * eyfgrid(j,i,n,m))
         eygrid(j,i) = eygrid(j,i)   & 
           + dble(ympole(n,m) * eyfgrid(j,i,n,m))
         ezgrid(j,i) = ezgrid(j,i)   & 
           + dble(zmpole(n,m) * eyfgrid(j,i,n,m))
         do m = 1,n
         exgrid(j,i) = exgrid(j,i)   & 
           + 2*dble(xmpole(n,m) * eyfgrid(j,i,n,m))
         eygrid(j,i) = eygrid(j,i)   & 
           + 2*dble(ympole(n,m) * eyfgrid(j,i,n,m))
         ezgrid(j,i) = ezgrid(j,i)   & 
           + 2*dble(zmpole(n,m) * eyfgrid(j,i,n,m))     
         enddo
        enddo  
       enddo
      enddo

      write(53,*) "istep =", istep
      write(53,*) "equal angle x y zgrid"  
      call log_xyz_grid(53, nsamp, nphi, exgrid, eygrid, ezgrid)
      mid = (nsamp+1)/2
      shortaxis = exgrid(1,mid)**2 + eygrid(1,mid)**2 + ezgrid(1,mid)**2
      shortaxis = dsqrt(shortaxis)
      longaxis = ezgrid(1,1)
      write(56,*) "istep =",istep,longaxis,shortaxis,longaxis/shortaxis
   
 end SUBROUTINE



 SUBROUTINE mod_eqsamp_eval_grid_xyz_pair(istep, &
  x1mpole, y1mpole, z1mpole, x2mpole, y2mpole, z2mpole)
  implicit real *8 (a-h,o-z)     
  integer :: istep
  complex *16 :: x1mpole(0:nterms, 0:nterms),y1mpole(0:nterms, 0:nterms),&
  z1mpole(0:nterms, 0:nterms)       
  complex *16 :: x2mpole(0:nterms, 0:nterms),y2mpole(0:nterms, 0:nterms),&
  z2mpole(0:nterms, 0:nterms)     

     do i = 1, nsamp
       do j = 1, nphi
        ex1grid(j,i) = 0d0
        ey1grid(j,i) = 0d0
        ez1grid(j,i) = 0d0
        do n = 0,nterms
         m = 0
         ex1grid(j,i) = ex1grid(j,i)   & 
           + dble(x1mpole(n,m) * eyfgrid(j,i,n,m))
         ey1grid(j,i) = ey1grid(j,i)   & 
           + dble(y1mpole(n,m) * eyfgrid(j,i,n,m))
         ez1grid(j,i) = ez1grid(j,i)   & 
           + dble(z1mpole(n,m) * eyfgrid(j,i,n,m))
         do m = 1,n
         ex1grid(j,i) = ex1grid(j,i)   & 
           + 2*dble(x1mpole(n,m) * eyfgrid(j,i,n,m))
         ey1grid(j,i) = ey1grid(j,i)   & 
           + 2*dble(y1mpole(n,m) * eyfgrid(j,i,n,m))
         ez1grid(j,i) = ez1grid(j,i)   & 
           + 2*dble(z1mpole(n,m) * eyfgrid(j,i,n,m))     
         enddo
        enddo  
       enddo
      enddo

      do i = 1, nsamp
        do j = 1, nphi
         ex2grid(j,i) = 0d0
         ey2grid(j,i) = 0d0
         ez2grid(j,i) = 0d0
         do n = 0,nterms
          m = 0
          ex2grid(j,i) = ex2grid(j,i)   & 
            + dble(x2mpole(n,m) * eyfgrid(j,i,n,m))
          ey2grid(j,i) = ey2grid(j,i)   & 
            + dble(y2mpole(n,m) * eyfgrid(j,i,n,m))
          ez2grid(j,i) = ez2grid(j,i)   & 
            + dble(z2mpole(n,m) * eyfgrid(j,i,n,m))
          do m = 1,n
          ex2grid(j,i) = ex2grid(j,i)   & 
            + 2*dble(x2mpole(n,m) * eyfgrid(j,i,n,m))
          ey2grid(j,i) = ey2grid(j,i)   & 
            + 2*dble(y2mpole(n,m) * eyfgrid(j,i,n,m))
          ez2grid(j,i) = ez2grid(j,i)   & 
            + 2*dble(z2mpole(n,m) * eyfgrid(j,i,n,m))     
          enddo
         enddo  
        enddo
       enddo

    write(153,*) "istep =", istep, "equal distributed theta: x y zgrid drop1"  
    call log_xyz_grid(153, nsamp, nphi, ex1grid, ey1grid, ez1grid)
    call log_xyz_grid(154, nsamp, nphi, ex1grid, ey1grid, ez1grid)
    write(253,*) "istep =", istep, "equal distributed theta: x y zgrid drop2"  
    call log_xyz_grid(253, nsamp, nphi, ex2grid, ey2grid, ez2grid)
    call log_xyz_grid(254, nsamp, nphi, ex2grid, ey2grid, ez2grid)

 end SUBROUTINE




 SUBROUTINE mod_eqsamp_eval_grid_fcomplex(mpole, fgrid)
  implicit real *8 (a-h,o-z)     
  
  complex *16 :: mpole(0:nterms, -nterms:nterms)
  complex *16 :: fgrid(nphi, nsamp)
  
     do i = 1, nsamp
       do j = 1, nphi
        fgrid(j,i) = 0d0
        do n = 0,nterms
         do m = -n,n
         fgrid(j,i) = fgrid(j,i) + mpole(n,m) * eyfgrid(j,i,n,m)
         enddo
        enddo  
       enddo
      enddo
   
 end SUBROUTINE


  SUBROUTINE mod_eqsamp_eval_grid_freal(mpole, fgrid)
  implicit real *8 (a-h,o-z)     
  
  complex *16 :: mpole(0:nterms, 0:nterms)
  real *8 :: fgrid(nphi, nsamp)
  
     do i = 1, nsamp
       do j = 1, nphi
        fgrid(j,i) = 0d0
        do n = 0,nterms
         m = 0
         fgrid(j,i) = fgrid(j,i) + dble(mpole(n,m) * eyfgrid(j,i,n,m))
         do m = 1,n
         fgrid(j,i) = fgrid(j,i) + 2*dble(mpole(n,m) * eyfgrid(j,i,n,m))
         enddo
        enddo  
       enddo
      enddo
   
 end SUBROUTINE  
 
end module mod_eqsamp
