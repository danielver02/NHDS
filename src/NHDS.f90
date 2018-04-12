! NHDS is a numerical solver to the linear hot-plasma dispersion relation
! Copyright (C) 2018 Daniel Verscharen (d.verscharen@ucl.ac.uk)
!
!This file is part of NHDS.
!
!    NHDS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    NHDS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with NHDS.  If not, see <http://www.gnu.org/licenses/>.

module input_params
  implicit none
  save
  integer :: numspec,numiter,nmax,mmax,ampl_mode
  integer :: vxsteps,vysteps,vzsteps,timesteps,num_periods
  double precision :: vtherm(10),alpha(10),Omega(10),ell(10),vdrift(10),mass(10),charge(10),beta(10),density(10)
  double precision :: vAc,det_D_threshold,Bessel_zero,ampl,theta
  double precision :: Bessel_zero_deltaf,vxrange(2),vyrange(2),vzrange(2)
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
x=0.009d0-0.0001d0*uniti

open(unit=10,file='output.dat',status='replace',action='write') 


do i=10,210,1
   kz=0.001d0*i
   kperp=kz*tan(theta*M_PI/180.d0)

  call newton_method(x,kz,kperp,quality)
  call calc_polarization(pol,polz,x,kz,kperp)
  call waveenergy(energy,x,kz,kperp,gamma_contribution)

  call calc_xi(xi,1,pol,polz,x,kz,kperp)	! second parameter is index of species 
  call calc_fluctRj(deltaRj,dUperpx,dUperpy,dUpar,1,pol,polz,x,kz,kperp) ! fifth parameter is index of species 

  write (*,*)  kz,real(x),aimag(x),real(pol),aimag(pol),real(polz),aimag(polz)!,energy,quality
  write (10,*) kz,real(x),aimag(x),real(pol),aimag(pol),real(polz),aimag(polz),energy,quality

! if (i.EQ.10) call write_delta_f(1,kz,x,pol,polz) ! first parameter is index of species   
enddo
 
close(10)



end program

