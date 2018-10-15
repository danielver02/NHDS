module input_params
  implicit none
  save
  integer :: numspec,numiter,nmax,mmax,ampl_mode
  integer :: vxsteps,vysteps,vzsteps,timesteps,num_periods,ksteps,theta_steps
  double precision :: vtherm(10),alpha(10),Omega(10),ell(10),vdrift(10),mass(10),charge(10),beta(10),density(10)
  double precision :: vAc,det_D_threshold,Bessel_zero,ampl,theta_range(2),krange(2)
  double precision :: Bessel_zero_deltaf,vxrange(2),vyrange(2),vzrange(2)
  double complex :: initial_guess
  double complex, parameter ::  uniti=(0.d0,1.d0)
  double precision,parameter :: M_PI=3.141592654d0
  logical :: output_warning,damping,periods,const_r,output_mom,output_EB,kth_file
  character*20 :: kth_filename
end module input_params


