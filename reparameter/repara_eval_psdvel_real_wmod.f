      SUBROUTINE eval_para_vel_real(nqterms, nqtheta, nqphi,
     $ idrop,   
     $ n_cutoff, nup, xtmpole, ytmpole, ztmpole,
     $ dthdtau, dphdtau)
      use OMP_LIB
      use mod_upsamp
ccccccccccccccccccccccccccccccccc
c      evaluate pseudo-velocity w and dthdtau, dphdtau
ccccccccccccccccccccccccccccccccc
      
       !variation of quality measure based on adaptive cut-off
       !with respect to grid positions
       
       !! evaluate n-cut-off
      implicit real *8 (a-h,o-z)
      !! input:
      integer :: nqterms, nqtheta, nqphi, n_cutoff, nup
      complex *16 :: xtmpole(0:nqterms, 0:nqterms),
     $ ytmpole(0:nqterms, 0:nqterms), ztmpole(0:nqterms, 0:nqterms) 
      !! output
      real* 8 :: dthdtau(nqphi, nqtheta),dphdtau(nqphi, nqtheta)


      real* 8, ALLOCATABLE :: dxdth(:,:,:), dxdph(:,:,:)   
      real* 8, ALLOCATABLE:: Geom_E(:,:), Geom_F(:,:), Geom_G(:,:)
      real* 8, ALLOCATABLE:: Geom_w(:,:), norm(:,:,:)
      real* 8, ALLOCATABLE :: var_quaMes(:,:,:), psd_vel(:,:,:)
      real* 8 :: psd_vel_mes
      DOUBLE PRECISION :: START, END  

      allocate(dxdth(nqphi,nqtheta,3), dxdph(nqphi,nqtheta,3))
      allocate(Geom_E(nqphi,nqtheta))
      allocate(Geom_F(nqphi,nqtheta))
      allocate(Geom_G(nqphi,nqtheta))
      allocate(Geom_W(nqphi,nqtheta), norm(nqphi,nqtheta,3))
      allocate(var_quaMes(nqphi,nqtheta,3), psd_vel(nqphi,nqtheta,3))
  
c$OMP PARALLEL 
c$OMP do schedule(static) private(i,n,m)  
      do j=1,nqphi
        do i=1,nqtheta
          dxdth(j,i,1)=0d0
          dxdth(j,i,2)=0d0
          dxdth(j,i,3)=0d0
          dxdph(j,i,1)=0d0
          dxdph(j,i,2)=0d0
          dxdph(j,i,3)=0d0
          do n = 0,nqterms
            m = 0
            dxdth(j,i,1)=dxdth(j,i,1)+
     $          dble(xtmpole(n,m)*upr_dydth(j,i,n,m))
            dxdth(j,i,2)=dxdth(j,i,2)+
     $          dble(ytmpole(n,m)*upr_dydth(j,i,n,m))
            dxdth(j,i,3)=dxdth(j,i,3)+
     $          dble(ztmpole(n,m)*upr_dydth(j,i,n,m))
            dxdph(j,i,1)=dxdph(j,i,1)+
     $          dble(xtmpole(n,m)*upr_dydph(j,i,n,m))
            dxdph(j,i,2)=dxdph(j,i,2)+
     $          dble(ytmpole(n,m)*upr_dydph(j,i,n,m))
            dxdph(j,i,3)=dxdph(j,i,3)+
     $          dble(ztmpole(n,m)*upr_dydph(j,i,n,m))
            do m = 1, n
cc        first order derivatives
               dxdth(j,i,1)=dxdth(j,i,1)+
     $            2*dble(xtmpole(n,m)*upr_dydth(j,i,n,m))
               dxdth(j,i,2)=dxdth(j,i,2)+
     $            2*dble(ytmpole(n,m)*upr_dydth(j,i,n,m))
               dxdth(j,i,3)=dxdth(j,i,3)+
     $            2*dble(ztmpole(n,m)*upr_dydth(j,i,n,m))
               dxdph(j,i,1)=dxdph(j,i,1)+
     $            2*dble(xtmpole(n,m)*upr_dydph(j,i,n,m))
               dxdph(j,i,2)=dxdph(j,i,2)+
     $            2*dble(ytmpole(n,m)*upr_dydph(j,i,n,m))
               dxdph(j,i,3)=dxdph(j,i,3)+
     $            2*dble(ztmpole(n,m)*upr_dydph(j,i,n,m))
            enddo
          enddo
 
          Geom_E(j,i)=dxdth(j,i,1)*dxdth(j,i,1)
     $    + dxdth(j,i,2)*dxdth(j,i,2) + dxdth(j,i,3)*dxdth(j,i,3)
          Geom_F(j,i)=dxdth(j,i,1)*dxdph(j,i,1)
     $    + dxdth(j,i,2)*dxdph(j,i,2) + dxdth(j,i,3)*dxdph(j,i,3)
          Geom_G(j,i)=dxdph(j,i,1)*dxdph(j,i,1)
     $    + dxdph(j,i,2)*dxdph(j,i,2) + dxdph(j,i,3)*dxdph(j,i,3) 
          Geom_W(j,i)=Geom_E(j,i)*Geom_G(j,i)-Geom_F(j,i)*Geom_F(j,i)
          Geom_W(j,i)=dsqrt(Geom_W(j,i))

        norm(j,i,1)=dxdth(j,i,2)*dxdph(j,i,3)-dxdth(j,i,3)*dxdph(j,i,2)
        norm(j,i,2)=dxdth(j,i,3)*dxdph(j,i,1)-dxdth(j,i,1)*dxdph(j,i,3)
        norm(j,i,3)=dxdth(j,i,1)*dxdph(j,i,2)-dxdth(j,i,2)*dxdph(j,i,1)
          norm(j,i,1) = norm(j,i,1) / Geom_W(j,i)
          norm(j,i,2) = norm(j,i,2) / Geom_W(j,i)
          norm(j,i,3) = norm(j,i,3) / Geom_W(j,i)             
       enddo
      enddo
c$OMP end do
c$OMP END PARALLEL    

      psd_vel_mes = 0d0 
c$OMP PARALLEL 
c$OMP do schedule(static) private(i,n,m,tmp,tmp_mes,k,b1,b2,
c$OMP& Ainver11,Ainver12,Ainver21,Ainver22, denomin)  
c$OMP& reduction(+:psd_vel_mes)
      do j=1,nqphi
        do i=1,nqtheta

          var_quaMes(j,i,1) = 0d0     
          var_quaMes(j,i,2) = 0d0   
          var_quaMes(j,i,3) = 0d0        
          do n = n_cutoff,nqterms  !! perfect low-pass filter
            m = 0
            var_quaMes(j,i,1) = var_quaMes(j,i,1) + 
     $        dble(xtmpole(n,m)*upr_sphfgrid(j,i,n,m))             
            var_quaMes(j,i,2) = var_quaMes(j,i,2) + 
     $        dble(ytmpole(n,m)*upr_sphfgrid(j,i,n,m))   
            var_quaMes(j,i,3) = var_quaMes(j,i,3) + 
     $        dble(ztmpole(n,m)*upr_sphfgrid(j,i,n,m))   
            do m = 1, n
              var_quaMes(j,i,1) = var_quaMes(j,i,1) + 
     $         2*dble(xtmpole(n,m)*upr_sphfgrid(j,i,n,m))             
              var_quaMes(j,i,2) = var_quaMes(j,i,2) + 
     $         2*dble(ytmpole(n,m)*upr_sphfgrid(j,i,n,m))   
              var_quaMes(j,i,3) = var_quaMes(j,i,3) + 
     $         2*dble(ztmpole(n,m)*upr_sphfgrid(j,i,n,m))             
            enddo
          enddo       
       
      !!scalar product
       tmp = 0d0
       do k = 1,3
         tmp = tmp + var_quaMes(j,i,k) * norm(j,i,k) 
       enddo
       tmp_mes = 0d0
       do k = 1,3
         psd_vel(j,i,k) =  -var_quaMes(j,i,k) + norm(j,i,k) * tmp
         tmp_mes = tmp_mes + psd_vel(j,i,k)**2
       enddo
       psd_vel_mes =  psd_vel_mes + dsqrt(tmp_mes) 

!!!debug: check the sign and mag of grid points near poles
       iwrite = 1000*idrop + 501
       if (i.eq.nqtheta) then
          write(iwrite,*) i,j,psd_vel(j,i,1),psd_vel(j,i,2),psd_vel(j,i,3)
       endif
!!!debug.

       !! scalar product
       b1 = 0d0
       b2 = 0d0
       Ainver11 = 0d0
       Ainver12 = 0d0
       Ainver22 = 0d0
       do k = 1, 3
         b1 = b1 + psd_vel(j,i,k) * dxdth(j,i,k)
         b2 = b2 + psd_vel(j,i,k) * dxdph(j,i,k)
         Ainver11 = Ainver11 + dxdph(j,i,k) * dxdph(j,i,k)
         Ainver12 = Ainver12 - dxdph(j,i,k) * dxdth(j,i,k)
         Ainver22 = Ainver22 + dxdth(j,i,k) * dxdth(j,i,k)
       enddo
       denomin = Ainver11*Ainver22 - Ainver12*Ainver12
       Ainver11 = Ainver11 / denomin
       Ainver12 = Ainver12 / denomin
       Ainver22 = Ainver22 / denomin
       Ainver21 = Ainver12
       
       dthdtau(j,i) = Ainver11 * b1 + Ainver12 * b2
       dphdtau(j,i) = Ainver21 * b1 + Ainver22 * b2
 !      write(48,*) "j,i,dthdtau,dphdtau",j,i,dthdtau(j,i),dphdtau(j,i)
      
       enddo
      enddo
c$OMP end do
c$OMP END PARALLEL  

      psd_vel_mes = psd_vel_mes / dble(nqtheta*nqphi)
      iwrite = idrop * 100 + 50
      write(iwrite,*) "psd_vel_mes = ", psd_vel_mes
      END SUBROUTINE
ccccccccccccccccccccccccccccccccc
