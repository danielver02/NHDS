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
integer :: j,i
double complex :: x 
double precision, allocatable, dimension(:) :: kk,theta
double complex, allocatable, dimension(:) :: guesses
double precision :: dk, dth
character*20 :: filename

dk=0.; dth=0.
call get_command_argument(1,filename); filename=adjustl(filename)
call set_parameters(filename)

allocate(kk(ksteps)) 
if (kth_file) then
   ! if reading k,th values from file
   ! ksteps=theta_steps=number of (k,th) pairs in input file
   allocate(theta(ksteps),guesses(ksteps))
   kth_filename=adjustl(kth_filename)
   open(unit=793,file=kth_filename)
   read(793,*) kk
   read(793,*) theta
   read(793,*) guesses
   close(793)
else
   allocate(theta(theta_steps))
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

open(unit=10,file='output_'//trim(filename)//'.dat',status='replace',action='write')
write (10,*) '# kk, theta, w, gamma, Re(iEy/Ex), Im(iEy/Ex), Re(Ez/Ex), Im(Ez/Ex), energy, quality'

if (output_mom) then
   open(unit=11,file='Plasma_'//trim(filename)//'.dat',status='replace',action='write')
   write (11,*) '# kk, theta, Re(dn), Im(dn), Re(dUperpx), Im(dUperpx), Re(dUperpy), Im(dUperpy), Re(dUparl), Im(dUparl)'
endif
if (output_EB) then
   open(unit=12,file='EB_'//trim(filename)//'.dat',status='replace',action='write')
   write (14,*) '# kk, theta, Re(Ex), Im(Ex), Re(Ey), Im(Ey), Re(Ez), Im(Ez), Re(Bx), Im(Bx), Re(By), Im(By), Re(Bz), Im(Bz)'
endif

if (kth_file) then
   do i=1,ksteps
      x = guesses(i) 
      call compute(kk(i),theta(i),x,output_mom,output_EB)
   enddo
else
   do j=1,theta_steps
      x=initial_guess ! guess w in units of wci:
      do i=1,ksteps
         call compute(kk(i),theta(j),x, output_mom,output_EB)
         if (i .eq. 1) initial_guess=x
      enddo
   enddo
endif
close(10); close(11); close(12)
end program

subroutine compute(kk,theta,x,outputm,outputeb)
   use input_params
   double precision, intent(in) :: kk, theta
   double complex, intent(inout) :: x
   logical, intent(in) :: outputm, outputeb
   double complex :: pol,polz,xi,deltaRj, Ek(3), Bk(3)
   double complex :: dUperpx,dUperpy,dUpar
   double precision :: kperp,kz,energy,gamma_contribution(10)
   real :: quality
   integer :: j
   kz=kk*cos(theta*M_PI/180.d0)
   kperp=kk*sin(theta*M_PI/180.d0)
   
   call newton_method(x,kz,kperp,quality)
   call calc_polarization(pol,polz,x,kz,kperp)
   call waveenergy(energy,x,kz,kperp,gamma_contribution)
   
!  if (i .eq. 1) initial_guess=x
   
   if (outputm) then
     do j=1,numspec
       call calc_xi(xi,j,pol,polz,x,kz,kperp)	
       ! second parameter is index of species
       call calc_fluctRj(deltaRj,dUperpx,dUperpy,dUpar,j,pol,polz,x,kz,kperp)
       ! fifth parameter is index of species
       write(11,'(10F14.5)') kk, theta, real(xi), aimag(xi), real(dUperpx), &
                            aimag(dUperpx), real(dUperpy), aimag(dUperpy), &
                            real(dUpar), aimag(dUpar)
     enddo
     write(11,*) ''
   endif
   
   if (outputeb) then
      Ek(1)=1.d0
      Ek(2)=pol*uniti
      Ek(3)=polz
      
      Bk(1)=-kz*Ek(2)/x/vAc
      Bk(2)=(kz*Ek(1)-kperp*Ek(3))/x/vAc
      Bk(3)=kperp*Ek(2)/x/vAc
      write(12,'(14F14.5)') kk, theta,   real(Ek(1)), aimag(Ek(1)), real(Ek(2))&
                          ,aimag(Ek(2)), real(Ek(3)), aimag(Ek(3)), real(Bk(1))&
                          ,aimag(Bk(1)), real(Bk(2)), aimag(Bk(2)), real(Bk(3))&
                          ,aimag(Bk(3))
   endif
   
!  write (*, '(10F9.5)') kk(i),theta(j),real(x),aimag(x),real(pol),aimag(pol),&
!                        real(polz),aimag(polz),energy,quality
   write (*, '(3F10.5,3E14.5)' ) kk, theta, real(x), aimag(x), energy, quality
   write (10,'(3F10.5,7E14.5)') kk, theta, real(x), aimag(x), real(pol), aimag(pol),&
                          real(polz), aimag(polz), energy, quality

! if (i.EQ.10) call write_delta_f(1,kz,x,pol,polz) ! first parameter is index of species
end subroutine compute
