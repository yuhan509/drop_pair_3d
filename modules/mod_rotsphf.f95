MODULE mod_rotsphf
 use mod_assign
 use mod_orgHelper
 
 complex *16, ALLOCATABLE :: rotmat(:,:,:,:), sphfgrid(:,:,:,:), & 
       rsphfgrids(:,:,:,:,:,:) 
 SAVE    
 
 CONTAINS
 
 !! must call mod_orgHelper_init first 
 SUBROUTINE mod_rotsphf_init()
     implicit real *8 (a-h,o-z)

      allocate(sphfgrid(nphi,ntheta,0:nterms,0:nterms))
      nrot = nphi
      nbeta = ntheta
      allocate(rotmat(0:nterms,-nterms:nterms,-nterms:nterms,nbeta))
      allocate(rsphfgrids(nphi,ntheta,nrot,nbeta,  &
       0:nterms,0:nterms))
     
      call rota_sphf(nterms, ntheta, nphi, nbeta, theta, nrot, &
      ctheta, whts, ynms, zwsave, &
      rotmat, sphfgrid, rsphfgrids) 

 END SUBROUTINE mod_rotsphf_init

END MODULE mod_rotsphf
