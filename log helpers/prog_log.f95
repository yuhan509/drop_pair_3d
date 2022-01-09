SUBROUTINE prog_log(iwrite,dt)
use mod_assign
!!! used for screen shot for interuppted prog
integer :: iwrite
real*8 :: dt
      write(iwrite,*) "All prameters"
      write(iwrite,*)  "nup_repara =", nup_repara, "ratio_cutoff =", ratio_cutoff,  & 
        "dtau =", dtau, "i_repara_max =", i_repara_max
      write(iwrite,*)  "nstep =", nstep,   & 
        "dt =", dt, "Ca =", Ca
      write(iwrite,*) "nterms =", nterms, "ntheta =", ntheta,  & 
        "maxup =", maxup, " maxup_crv =", maxup_crv
      write(iwrite,*) "initial pair center_dist =", center_dist
      write(iwrite,*) "initial time step = ",dt
      
      write(iwrite,*) "electric field angle with z axis =", alpha
      write(iwrite,*) "nintpol=",nintpol, " hintpol=", hintpol," h_far=",h_far
      !write(iwrite,*) "tol_squ_mat =", 1d-25!, "tol_squ_rep =", 2d-28
      write(iwrite,*) "march_tol = ", march_tol, "tmax=", tmax
  

RETURN
END SUBROUTINE
