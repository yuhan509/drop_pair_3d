module mod_assign
 
 integer,PARAMETER :: nterms = 19
 integer,PARAMETER :: maxup = 7, maxup_crv = 4, nup_repara = 2
 integer,PARAMETER :: mterms = maxup_crv * nterms, &
  mxterms = maxup * nterms, nrterms = nup_repara * nterms
 integer :: mtheta, mphi, mxtheta, mxphi, nrtheta, nrphi
 integer :: ntheta, nphi 
  !! parameters control simulation
 integer :: nintpol, nstep
 real *8 :: Ca, alpha ,hintpol, h_far, march_tol, tmax
 real *8 :: center_dist
 integer :: i_repara_max
 real *8 :: ratio_cutoff, dtau
 save
 CONTAINS
 
end module mod_assign
