module mod_upsamp
 use mod_assign
 
 real *8, ALLOCATABLE,dimension(:,:) :: up_test
!!cccccccccc pre-calculated arrays for upsampling ccccccccccc
 real* 8, ALLOCATABLE, dimension(:,:) :: up_ctheta, up_whts, up_theta
 real* 8, ALLOCATABLE :: up_ynms(:,:,:,:) 
 complex *16, ALLOCATABLE, dimension(:,:) :: up_zwsave, up_dwsave
 complex *16, ALLOCATABLE, dimension(:,:,:,:,:) :: up_sphfgrid, &
 up_dydph, up_dydth, up_d2ydth2, up_d2ydph2, up_d2ydthdph
 complex *16, ALLOCATABLE, dimension(:,:,:,:) :: upr_sphfgrid, &
 upr_dydph, upr_dydth
!!ccccccc pre-calculate for truncate mpole remedy in upsampling cccc    
 complex *16, ALLOCATABLE, dimension(:,:,:,:) :: up_yfngrid, yfngrid, &
 upx_sphfgrid

 CONTAINS 
 
 SUBROUTINE init_test
  implicit real *8 (a-h,o-z)
  allocate(up_test(nterms, 2 * nterms))
 
  do i = 1, nterms
    do j = 1, 2 * nterms
      up_test(i,j) = i + j
    enddo
  enddo
  END SUBROUTINE
  
  subroutine mod_upsamp_mterms_init()
  
    call grid_dim_init(mterms, mtheta, mphi)
    allocate(up_sphfgrid(mphi,mtheta,0:nterms,0:nterms,maxup_crv))
    allocate(up_dydph(mphi,mtheta,0:nterms,0:nterms,maxup_crv))
    allocate(up_dydth(mphi,mtheta,0:nterms,0:nterms,maxup_crv))
    allocate(up_d2ydth2(mphi,mtheta,0:nterms,0:nterms,maxup_crv))      
    allocate(up_d2ydph2(mphi,mtheta,0:nterms,0:nterms,maxup_crv))  
    allocate(up_d2ydthdph(mphi,mtheta,0:nterms,0:nterms,maxup_crv))
    
    call grid_dim_init(mxterms, mxtheta, mxphi)
    allocate(upx_sphfgrid(mxphi,mxtheta,0:nterms,0:nterms))

    allocate(up_ctheta(mxtheta,maxup), up_whts(mxtheta,maxup), up_theta(mxtheta,maxup))
    allocate(up_ynms(0:mxterms,0:mxterms,mxtheta/2+1,maxup))
    allocate(up_zwsave(4*mxphi+15,maxup),up_dwsave(4*mxphi+15,maxup))
    
    call grid_dim_init(nrterms, nrtheta, nrphi)    
    allocate(upr_sphfgrid(nrphi,nrtheta,0:nrterms,0:nrterms))
    allocate(upr_dydph(nrphi,nrtheta,0:nrterms,0:nrterms))
    allocate(upr_dydth(nrphi,nrtheta,0:nrterms,0:nrterms))
    

      !! all the pre-calculated arrays for upsampling (curvatrue evaluation)
      !! ctheta, whts, ynms, dwsave, zwsave(for sphfgrid)
      !! sphfgrid, dydph, dydth, d2ydth2, d2ydph2, d2ydthdph
    call upsamp_crv_pre_cal()
    call upsamp_upx_pre_cal()
    call upsamp_repara_pre_cal()
   
  end SUBROUTINE mod_upsamp_mterms_init

  !! called after init ntheta, nphi
  SUBROUTINE mod_upsamp_org_grid_init(theta)
     implicit real *8 (a-h,o-z)     
     real*8 :: theta(ntheta)
     real *8, ALLOCATABLE :: phigrid(:,:), thetagrid(:,:)

     allocate(phigrid(nphi,ntheta),thetagrid(nphi,ntheta))
     PI = 4.D0*DATAN(1.D0)
      phi = 2*PI/nphi
      do i=1,ntheta
        do j=1,nphi
          phigrid(j, i) = (j-1)*phi
          thetagrid(j, i) = theta(i)
        enddo
      enddo    
    
    allocate(up_yfngrid(nphi,ntheta,0:mterms,-mterms:mterms))
    call sph_grid(mterms,nphi,ntheta,thetagrid,phigrid,up_yfngrid)  
    deallocate(phigrid,thetagrid)
    
  end SUBROUTINE
  
end module mod_upsamp
