! This file is part of NHDS
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

subroutine waveenergy(energy,x,kz,kperp,gamma_contribution)
use input_params
implicit none
double complex :: ep(3,3),chi(3,3),chiplus(3,3),Ek(3),Bk(3),epprime(3,3), chi_a(10,3,3)
double complex :: Ekmult(3), Ekmultprime(3), eph(3,3), epprimeh(3,3)
double precision :: kz,kperp,energy,gamma_contribution(10)
double complex :: x,pol,polz,xr,dx
integer :: i,k,j

dx=1.d-7


do i=1,3
 do k=1,3
  ep(i,k)=0.d0
  epprime(i,k)=0.d0
 enddo
enddo

xr=x-uniti*aimag(x)


! evaluate at omega = omega_{kr}
do j=1,numspec
  call calc_chi(chiplus,j,kz,kperp,xr+dx)
  call calc_chi(chi,j,kz,kperp,xr-dx)
  do i=1,3
    do k=1,3
      ep(i,k)=ep(i,k)+chi(i,k)
      epprime(i,k)=epprime(i,k)+chiplus(i,k)
    enddo
  enddo
enddo


! derivative of ep at omega = omega_{kr}:
do i=1,3
 do k=1,3
  epprime(i,k)=(epprime(i,k)-ep(i,k))/(2.d0*dx)
  ep(i,k)=0.d0
 enddo
enddo


! evaluate ep at omega = omega_{kr}
do j=1,numspec
  call calc_chi(chi,j,kz,kperp,xr)
  do i=1,3
    do k=1,3
      !anti-hermitian part of the susceptibility:
  	  chi_a(j,i,k)=(chi(i,k)-(real(chi(k,i))-uniti*aimag(chi(k,i))))/(2.d0*uniti)
      ep(i,k)=ep(i,k)+chi(i,k)
    enddo
  enddo
enddo


!vAc is the Alfven speed in units of speed of light
do i=1,3
  ep(i,i)=ep(i,i)+vAc*vAc
enddo

call calc_polarization(pol,polz,x,kz,kperp)

! use all field components in units of E_kx:
Ek(1)=1.d0
Ek(2)=pol*uniti
Ek(3)=polz


Bk(1)=-kz*Ek(2)/x
Bk(2)=(kz*Ek(1)-kperp*Ek(3))/x
Bk(3)=kperp*Ek(2)/x



! calculate Hermitian part only:
do i=1,3
	do k=1,3
		eph(i,k)=(ep(i,k)+1.d0*real(ep(k,i))-uniti*aimag(ep(k,i)))/2.d0
		epprimeh(i,k)=(epprime(i,k)+1.d0*real(epprime(k,i))-uniti*aimag(epprime(k,i)))/2.d0
	enddo
enddo





! Now we have the dielectric tensor and the polarization. Let's calculate the wave energy

do i=1,3
	Ekmult(i)=0.d0
	Ekmultprime(i)=0.d0
  	do k=1,3
    	Ekmult(i)=Ekmult(i)+eph(i,k)*Ek(k)
    	Ekmultprime(i)=Ekmultprime(i)+epprimeh(i,k)*Ek(k)
	enddo
enddo


energy=0.d0
do i=1,3
	energy=energy+1.d0*real((1.d0*real(Ek(i))-uniti*aimag(Ek(i)))*Ekmult(i))
	energy=energy+xr*real((1.d0*real(Ek(i))-uniti*aimag(Ek(i)))*Ekmultprime(i))
	energy=energy+1.d0*real(Bk(i))**2+aimag(Bk(i))**2
enddo

energy=energy/(16.d0*M_PI*vAc*vAc)


! Determine the contribution to gamma from all species (see Eq. (4) in Quataert (1998))
gamma_contribution=0.d0

do j=1,numspec
	! chi_a and Ek are already correct from above.
	do k=1,3
		do i=1,3
			gamma_contribution(j) = gamma_contribution(j)+ (real(Ek(k))-uniti*aimag(Ek(k)))*chi_a(j,k,i)*Ek(i)
		enddo
	enddo
	gamma_contribution(j) = -xr*gamma_contribution(j) /(4.d0*vAc*vAc*energy*4.d0*M_PI)
enddo
end subroutine
