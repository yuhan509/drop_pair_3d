      subroutine shf_infi_area(iterms, itheta, iphi,
     $ ximpole, yimpole, zimpole, dydph, dydth,
     $ wigrid, nigrid)
     
      use mod_upsamp

      implicit real *8 (a-h,o-z)
      integer :: iterms, itheta, iphi, iup
      complex *16 :: ximpole(0:iterms,0:iterms),
     $ yimpole(0:iterms,0:iterms),zimpole(0:iterms,0:iterms)
      real *8 :: wigrid(iphi, itheta),
     $ nigrid(iphi, itheta,3)

      complex *16 ::  !! change to real *8
     $ dydph(iphi,itheta,0:iterms,0:iterms), 
     $ dydth(iphi,itheta,0:iterms,0:iterms)
      real* 8, ALLOCATABLE, dimension(:,:,:):: dxdth, dxdph
      real* 8, ALLOCATABLE:: Geom_E(:,:),Geom_F(:,:),Geom_G(:,:)

      allocate(dxdth(iphi,itheta,3))
      allocate(dxdph(iphi,itheta,3))

      allocate(Geom_E(iphi,itheta))
      allocate(Geom_F(iphi,itheta))
      allocate(Geom_G(iphi,itheta))

      do i=1,itheta
        do j=1,iphi
          do k = 1,3
cc        first order derivatives              
          dxdth(j,i,k) = 0d0                
          dxdph(j,i,k) = 0d0        
          enddo
   
          do l = 0,iterms
            m = 0
cc        first order derivatives
        dxdth(j,i,1)=dxdth(j,i,1)+dble(ximpole(l,m)*dydth(j,i,l,m))
        dxdth(j,i,2)=dxdth(j,i,2)+dble(yimpole(l,m)*dydth(j,i,l,m))
        dxdth(j,i,3)=dxdth(j,i,3)+dble(zimpole(l,m)*dydth(j,i,l,m))
          
        dxdph(j,i,1)=dxdph(j,i,1)+dble(ximpole(l,m)*dydph(j,i,l,m))
        dxdph(j,i,2)=dxdph(j,i,2)+dble(yimpole(l,m)*dydph(j,i,l,m))
        dxdph(j,i,3)=dxdph(j,i,3)+dble(zimpole(l,m)*dydph(j,i,l,m))         
            do m = 1, l 
cc        first order derivatives
        dxdth(j,i,1)=dxdth(j,i,1)+ 2*dble(ximpole(l,m)*dydth(j,i,l,m))
        dxdth(j,i,2)=dxdth(j,i,2)+ 2*dble(yimpole(l,m)*dydth(j,i,l,m))
        dxdth(j,i,3)=dxdth(j,i,3)+ 2*dble(zimpole(l,m)*dydth(j,i,l,m))

        dxdph(j,i,1)=dxdph(j,i,1)+ 2*dble(ximpole(l,m)*dydph(j,i,l,m))
        dxdph(j,i,2)=dxdph(j,i,2)+ 2*dble(yimpole(l,m)*dydph(j,i,l,m))
        dxdph(j,i,3)=dxdph(j,i,3)+ 2*dble(zimpole(l,m)*dydph(j,i,l,m))

            enddo
          enddo

          Geom_E(j,i)=dxdth(j,i,1)*dxdth(j,i,1)
     $    + dxdth(j,i,2)*dxdth(j,i,2) + dxdth(j,i,3)*dxdth(j,i,3)
          Geom_F(j,i)=dxdth(j,i,1)*dxdph(j,i,1)
     $    + dxdth(j,i,2)*dxdph(j,i,2) + dxdth(j,i,3)*dxdph(j,i,3)
          Geom_G(j,i)=dxdph(j,i,1)*dxdph(j,i,1)
     $    + dxdph(j,i,2)*dxdph(j,i,2) + dxdph(j,i,3)*dxdph(j,i,3) 
          wigrid(j,i)=Geom_E(j,i)*Geom_G(j,i) - Geom_F(j,i)**2
          wigrid(j,i)=dsqrt(wigrid(j,i))

      nigrid(j,i,1)=dxdth(j,i,2)*dxdph(j,i,3)-dxdth(j,i,3)*dxdph(j,i,2)
      nigrid(j,i,2)=dxdth(j,i,3)*dxdph(j,i,1)-dxdth(j,i,1)*dxdph(j,i,3)
      nigrid(j,i,3)=dxdth(j,i,1)*dxdph(j,i,2)-dxdth(j,i,2)*dxdph(j,i,1)
        nigrid(j,i,1) =nigrid(j,i,1) / wigrid(j,i)
        nigrid(j,i,2) =nigrid(j,i,2) / wigrid(j,i)
        nigrid(j,i,3) =nigrid(j,i,3) / wigrid(j,i)
        enddo
      enddo 

      END SUBROUTINE