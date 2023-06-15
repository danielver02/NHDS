! NHDS is a numerical solver to the linear hot-plasma dispersion relation
! Copyright (C) 2020 Daniel Verscharen (d.verscharen@ucl.ac.uk)
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
  integer :: numspec,numiter,nmax,mmax,ampl_mode,species_df
  integer :: vxsteps,vysteps,vzsteps,timesteps,num_periods,ksteps,theta_steps
  double precision :: vtherm(10),alpha(10),Omega(10),ell(10),vdrift(10),mass(10),charge(10),beta(10),density(10)
  double precision :: vAc,det_D_threshold,Bessel_zero,ampl,theta_range(2),krange(2)
  double precision :: Bessel_zero_deltaf,vxrange(2),vyrange(2),vzrange(2)
  double complex :: initial_guess
  double complex, parameter ::  uniti=(0.d0,1.d0)
  double precision,parameter :: M_PI=3.141592654d0
  logical :: output_warning,damping,periods,const_r,output_mom,output_EB,kth_file,output_df
  character*300 :: kth_filename
end module input_params

program NHDS
use input_params
implicit none
integer :: j,i,l
double complex :: x
double precision, allocatable, dimension(:) :: kk,theta
double precision :: dk, dth
character*300 :: filename

dk=0.; dth=0.
call get_command_argument(1,filename); filename=adjustl(filename)
call set_parameters(filename)

allocate(kk(ksteps))
if (kth_file) then
   ! if reading k,th values from file
   ! ksteps=theta_steps=number of (k,th) pairs in input file
   allocate(theta(ksteps))
   kth_filename=adjustl(kth_filename)
   open(unit=793,file=kth_filename)
   do i=1,ksteps
      read(793,*) kk(i), theta(i)
    enddo
   close(793)
else
   allocate(theta(theta_steps))
   ! If using range of k and theta, create the kk, theta arrays to use
   if (ksteps .gt. 1) dk = (krange(2)-krange(1))/(1.d0*ksteps-1.d0)
   if (theta_steps .gt. 1) dth = (theta_range(2)-theta_range(1))/(1.d0*theta_steps-1.d0)
   do i=1,ksteps
      kk(i)=krange(1)+dk*(1.d0*i-1.d0)
    enddo
   do i=1,theta_steps
      theta(i)=theta_range(1)+dth*(1.d0*i-1.d0)
   enddo
endif

do j=1,numspec
  Omega(j)=charge(j)/mass(j)
  ell(j)=sqrt(mass(j)/(density(j)*charge(j)*charge(j)))
  vtherm(j)=sqrt(beta(j)/(density(j)*mass(j)))
enddo

open(unit=10,file='output_'//trim(filename)//'_wave.dat',status='replace',action='write')
write (10,*) '# kk, theta, omega, gamma, Re(Ey/Ex), Im(Ey/Ex), Re(Ez/Ex), Im(Ez/Ex), energy, quality'

if (output_mom) then
  open(unit=11,file='output_'//trim(filename)//'_plasma.dat',status='replace',action='write')
   write (11,*) '# kk, theta, j, Re(dn), Im(dn), Re(dUperpx), Im(dUperpx), Re(dUperpy), Im(dUperpy), Re(dUpar), &
            Im(dUpar), Re(dpper), Im(dpperp), Re(dppar), Im(dppar), gamma_contribution, heating rate'
endif
if (output_EB) then
  open(unit=12,file='output_'//trim(filename)//'_EB.dat',status='replace',action='write')
   write (12,*) '# kk, theta, Re(Ex), Im(Ex), Re(Ey), Im(Ey), Re(Ez), Im(Ez), Re(Bx), Im(Bx), Re(By), Im(By), Re(Bz), Im(Bz)'
endif

if (kth_file) then
   x=initial_guess ! guess w in units of wci:
   do i=1,ksteps
      call compute(kk(i),theta(i),x,output_mom,output_EB)
   enddo
else
   do l=1,theta_steps
      x=initial_guess ! guess w in units of wci:
      do i=1,ksteps
         call compute(kk(i),theta(l),x, output_mom,output_EB)
         if (i .eq. 1) initial_guess=x
      enddo
   enddo
endif

close(10)
close(11)
close(12)

end program

subroutine compute(kk,theta,x,outputm,outputeb)
   use input_params
   double precision, intent(in) :: kk, theta
   double complex, intent(inout) :: x
   logical, intent(in) :: outputm, outputeb
   double complex :: Avec(3),xi,Ek(3), Bk(3)
   double complex :: dUperpx,dUperpy,dUpar,dpperp,dppar,Avec(3)
   double precision :: kperp,kz,energy,gamma_contribution(10), heating(10)
   double precision :: quality
   integer :: j
   kz=kk*cos(theta*M_PI/180.d0)
   kperp=kk*sin(theta*M_PI/180.d0)

   call newton_method(x,kz,kperp,quality)
   call calc_polarization(Avec,Ek,Bk,x,kz,kperp)
   call waveenergy(energy,x,kz,kperp,gamma_contribution,heating)

!  if (i .eq. 1) initial_guess=x

   if (outputm) then
     do j=1,numspec
       call calc_xi(xi,j,Avec,x,kz,kperp)
       ! second parameter is index of species
       call calc_fluctRj(dpperp,dppar,dUperpx,dUperpy,dUpar,j,Avec,x,kz,kperp)
       ! fifth parameter is index of species
       write(11,'(ES25.15E3,F15.10,I3,14ES25.15E3)') kk, theta, j, real(xi), aimag(xi), real(dUperpx), &
                            aimag(dUperpx), real(dUperpy), aimag(dUperpy), &
                            real(dUpar), aimag(dUpar), real(dpperp), aimag(dpperp), &
                            real(dppar), aimag(dppar), gamma_contribution(j), heating(j)
     enddo
   endif

   if (outputeb) then

      write(12,'(ES25.15E3,F15.10,12ES25.15E3)') kk, theta,   real(Ek(1)), aimag(Ek(1)), real(Ek(2))&
                          ,aimag(Ek(2)), real(Ek(3)), aimag(Ek(3)), real(Bk(1))&
                          ,aimag(Bk(1)), real(Bk(2)), aimag(Bk(2)), real(Bk(3))&
                          ,aimag(Bk(3))
   endif

!  write (*, '(10F9.5)') kk(i),theta(j),real(x),aimag(x),real(pol),aimag(pol),&
!                        real(polz),aimag(polz),energy,quality
   write (*, '(ES18.6E3,F10.5,3ES18.6E3)' ) kk, theta, real(x), aimag(x), quality
   write (10,'(ES25.15E3,F15.10,8ES25.15E3)') kk, theta, real(x), aimag(x), real(Ek(2)/EK(1)), aimag(Ek(2)/EK(1)),&
                          real(EK(3)/EK(1)), aimag(EK(3)/EK(1)), energy, quality


   if (output_df) call write_delta_f(species_df,kk,theta,x,Avec) ! first parameter is index of species


end subroutine compute
