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

program NHDS
use input_params
implicit none
integer :: j,i,k
double complex :: x,pol,polz,xi,deltaRj
double complex :: dUperpx,dUperpy,dUpar
double precision, allocatable, dimension(:) :: kk,theta
double precision :: dk, dth
double precision :: kperp,kz,energy,gamma_contribution(10)
real :: quality

dk=0.; dth=0.
call set_parameters()

allocate(kk(ksteps), theta(theta_steps))
if (kth_file) then
   ! if reading k,th values from file
   ! ksteps=theta_steps=number of (k,th) pairs in input file
   kth_filename=adjustl(kth_filename)
   !write(*,*) kth_filename
   open(unit=793,file=kth_filename)
   read(793,*) kk
   read(793,*) theta
   close(793)
else
   ! If using range of k and theta, create the kk, theta arrays to use
   if (ksteps .gt. 1) dk = (krange(2)-krange(1))/(1.d0*ksteps-1.d0)
   if (theta_steps .gt. 1) dth = (theta_range(2)-theta_range(1))/(1.d0*theta_steps-1.d0)
   do i=1,ksteps; kk(i)=krange(1)+dk*(1.d0*i-1.d0); enddo
   do i=1,theta_steps; theta(i)=theta_range(1)+dth*(1.d0*i-1.d0); enddo
endif

do j=1,numspec
  Omega(j)=charge(j)/mass(j)
  ell(j)=sqrt(mass(j)/(density(j)*charge(j)*charge(j)))
  vtherm(j)=sqrt(beta(j)/(density(j)*mass(j)))
enddo

open(unit=10,file='output.dat',status='replace',action='write')
write (10,*) '# kk, theta, w, gamma, Re(iEy/Ex), Im(iEy/Ex), Re(iEz/Ex), Im(iEz/Ex), energy, quality'

if (output_mom) then
   open(unit=11,file='Plasma.dat',status='replace',action='write')
   write (11,*) '# kk, theta, Re(dn), Im(dn), Re(dUperpx), Im(dUperpx), Re(dUperpy), Im(dUperpy), Re(dUparl), Im(dUparl)'
endif

do j=1,theta_steps

   ! Initial guess for frequency in units of proton gyrofrequency:
   x=initial_guess
   do i=1,ksteps

      kz=kk(i)*cos(theta(j)*M_PI/180.d0)
      kperp=kk(i)*sin(theta(j)*M_PI/180.d0)
      
      call newton_method(x,kz,kperp,quality)
      call calc_polarization(pol,polz,x,kz,kperp)
      call waveenergy(energy,x,kz,kperp,gamma_contribution)
      
      if (i .eq. 1) initial_guess=x
      
      if (output_mom) then
         do k=1,numspec
            call calc_xi(xi,k,pol,polz,x,kz,kperp)	! second parameter is index of species
            call calc_fluctRj(deltaRj,dUperpx,dUperpy,dUpar,k,pol,polz,x,kz,kperp) ! fifth parameter is index of species
            write(11,'(10F9.5)') kk(i), theta(j), real(xi), aimag(xi), real(dUperpx), aimag(dUperpx), real(dUperpy), aimag(dUperpy), real(dUpar), aimag(dUpar)
         enddo
         write(11,*) ''
      endif
      
      write (*, '(10F9.5)') kk(i),theta(j),real(x),aimag(x),real(pol),aimag(pol),real(polz),aimag(polz),energy,quality
      write (10,'(10F9.5)') kk(i),theta(j),real(x),aimag(x),real(pol),aimag(pol),real(polz),aimag(polz),energy,quality

! if (i.EQ.10) call write_delta_f(1,kz,x,pol,polz) ! first parameter is index of species
   enddo
enddo

close(10); close(11)
end program
