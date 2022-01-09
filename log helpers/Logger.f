      SUBROUTINE log_xyz_expansion_energy(iwrite,nterms,xmpole,ympole,
     $ zmpole)
      implicit real *8 (a-h,o-z)
      integer :: iwrite, nterms
      complex *16 :: xmpole(0:nterms,-nterms:nterms),
     $ ympole(0:nterms,-nterms:nterms), zmpole(0:nterms,-nterms:nterms)
      
      write(iwrite,*) "log_xyz_expansion_energy, iwrite =", iwrite
      write(iwrite,*) "nterms =", nterms
      write(iwrite,*) "n, En, Enx, Eny, Enz"

      do n = 0, nterms
        En = 0d0
        Enx = 0d0
        Eny = 0d0
        Enz = 0d0
        do m = -n, n           
          Enx = Enx + dsqrt(dble(xmpole(n,m)*dconjg(xmpole(n,m)))) 
          Eny = Eny + dsqrt(dble(ympole(n,m)*dconjg(ympole(n,m))))  
          Enz = Enz + dsqrt(dble(zmpole(n,m)*dconjg(zmpole(n,m))))
        enddo
        En = Enx + Eny + Enz
        write(iwrite,*) n, En, Enx, Eny, Enz
      enddo           
      END SUBROUTINE
c
c
c
      SUBROUTINE log_xyz_real_expan_energy(iwrite,nterms,xmpole,ympole,
     $ zmpole)
      implicit real *8 (a-h,o-z)
      integer :: iwrite, nterms
      complex *16 :: xmpole(0:nterms,0:nterms),
     $ ympole(0:nterms,0:nterms), zmpole(0:nterms,0:nterms)
      
      write(iwrite,*) "log_xyz_expansion_energy, iwrite =", iwrite
      write(iwrite,*) "nterms =", nterms
      write(iwrite,*) "n, En, Enx, Eny, Enz"

      do n = 0, nterms
         En = 0d0
        Enx = 0d0
        Eny = 0d0
        Enz = 0d0   
         m = 0
          Enx = Enx + dsqrt(dble(xmpole(n,m)*dconjg(xmpole(n,m)))) 
          Eny = Eny + dsqrt(dble(ympole(n,m)*dconjg(ympole(n,m))))  
          Enz = Enz + dsqrt(dble(zmpole(n,m)*dconjg(zmpole(n,m))))
        do m = 1, n           
          Enx = Enx + 2*dsqrt(dble(xmpole(n,m)*dconjg(xmpole(n,m)))) 
          Eny = Eny + 2*dsqrt(dble(ympole(n,m)*dconjg(ympole(n,m))))  
          Enz = Enz + 2*dsqrt(dble(zmpole(n,m)*dconjg(zmpole(n,m))))
        enddo
        En = Enx + Eny + Enz
        write(iwrite,*) n, En, Enx, Eny, Enz
      enddo           
      END SUBROUTINE
c
c
c
      SUBROUTINE log_cxyz_grid(iwrite,ntheta,nphi,xgrid,ygrid,zgrid)
      implicit real *8 (a-h,o-z)
      integer :: iwrite,ntheta,nphi
      complex *16 :: xgrid(nphi,ntheta),ygrid(nphi,ntheta),
     $ zgrid(nphi,ntheta)
      
      write(iwrite,*) "log_cxyz_grid, iwrite =", iwrite
      write(iwrite,*) "ntheta =",ntheta,"nphi =",nphi
      write(iwrite,*) "x,y,zgrid, & a value should be 1.0"

      do i = 1, ntheta
        do j = 1, nphi           
            write(iwrite,*) i,j,dble(xgrid(j,i)),dble(ygrid(j,i)),
     $ dble(zgrid(j,i)),dble(xgrid(j,i))**2 + dble(ygrid(j,i))**2 
     $ + dble(zgrid(j,i))**2/4   
        enddo
      enddo 
          
      END SUBROUTINE
      
      
      SUBROUTINE log_xyz_grid(iwrite,ntheta,nphi,xgrid,ygrid,zgrid)
      implicit real *8 (a-h,o-z)
      integer :: iwrite,ntheta,nphi
      real *8 :: xgrid(nphi,ntheta),ygrid(nphi,ntheta),
     $ zgrid(nphi,ntheta)
      
    !  write(iwrite,*) "log_xyz_grid, iwrite =", iwrite
    !  write(iwrite,*) "ntheta =",ntheta,"nphi =",nphi
    !  write(iwrite,*) "x,y,zgrid, & a value should be 1.0"

      do i = 1, ntheta
        do j = 1, nphi  
          rad = xgrid(j,i)**2 + ygrid(j,i)**2 + zgrid(j,i)**2
          rad = dsqrt(rad)           
          write(iwrite,*) i,j,xgrid(j,i),ygrid(j,i),zgrid(j,i)
        enddo
      enddo 
          
      END SUBROUTINE
