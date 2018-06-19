! NHDS is a numerical solver to the linear hot-plasma dispersion relation
! Copyright (C) 2018 Daniel Verscharen (d.verscharen@ucl.ac.uk)
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!The views and conclusions contained in the software and documentation are those
!of the authors and should not be interpreted as representing official policies,
!either expressed or implied, of the NHDS project.

module input_params
  implicit none
  save
  integer :: numspec,numiter,nmax,mmax,ampl_mode
  integer :: vxsteps,vysteps,vzsteps,timesteps,num_periods,kzsteps
  double precision :: vtherm(10),alpha(10),Omega(10),ell(10),vdrift(10),mass(10),charge(10),beta(10),density(10)
  double precision :: vAc,det_D_threshold,Bessel_zero,ampl,theta,kzrange(2)
  double precision :: Bessel_zero_deltaf,vxrange(2),vyrange(2),vzrange(2)
  double complex :: initial_guess
  double complex, parameter ::  uniti=(0.d0,1.d0)
  double precision,parameter :: M_PI=3.141592654d0
  logical :: output_warning,damping,periods,const_r
end module input_params

program NHDS
use input_params
implicit none
integer :: j,i
double complex :: x,pol,polz,xi,deltaRj
double complex :: dUperpx,dUperpy,dUpar
double precision :: kperp,kz,energy,gamma_contribution(10)
real :: quality


call set_parameters()

do j=1,numspec
  Omega(j)=charge(j)/mass(j)
  ell(j)=sqrt(mass(j)/(density(j)*charge(j)*charge(j)))
  vtherm(j)=sqrt(beta(j)/(density(j)*mass(j)))
enddo

! Initial guess for frequency in units of proton gyrofrequency:
x=initial_guess

open(unit=10,file='output.dat',status='replace',action='write')


do i=1,kzsteps

   kz=kzrange(1)+(kzrange(2)-kzrange(1))*(1.d0*i-1.d0)/(1.d0*kzsteps-1.d0)
   kperp=kz*tan(theta*M_PI/180.d0)

  call newton_method(x,kz,kperp,quality)
  call calc_polarization(pol,polz,x,kz,kperp)
  call waveenergy(energy,x,kz,kperp,gamma_contribution)

  call calc_xi(xi,1,pol,polz,x,kz,kperp)	! second parameter is index of species
  call calc_fluctRj(deltaRj,dUperpx,dUperpy,dUpar,1,pol,polz,x,kz,kperp) ! fifth parameter is index of species

  write (*,*)  kz,real(x),aimag(x),real(pol),aimag(pol),real(polz),aimag(polz),energy,quality
  write (10,*) kz,real(x),aimag(x),real(pol),aimag(pol),real(polz),aimag(polz),energy,quality

! if (i.EQ.10) call write_delta_f(1,kz,x,pol,polz) ! first parameter is index of species
enddo

close(10)



end program
