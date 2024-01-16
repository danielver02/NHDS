! This file is part of NHDS
! Copyright (C) 2024 Daniel Verscharen (d.verscharen@ucl.ac.uk)
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

subroutine waveenergy(energy,x,kz,kperp,gamma_contribution,heating)
use input_params
implicit none
double complex :: ep(3,3),chi(3,3),chiplus(3,3),Avec(3),Ek(3),Bk(3),epprime(3,3)
double complex, allocatable, dimension(:,:,:) :: chi_a
double complex :: Ekmult(3), Ekmultprime(3), eph(3,3), epprimeh(3,3)
double precision :: kz,kperp,energy,gamma_contribution(numspec), heating(numspec)
double complex :: x,xr,dx
integer :: i,k,j

allocate(chi_a(numspec,3,3))

xr=real(x)
dx=1.d-9*xr

ep=0.d0
epprime=0.d0

! evaluate at omega = omega_{kr}
do j=1,numspec
  call calc_chi(chiplus,j,kz,kperp,xr+dx)
  call calc_chi(chi,j,kz,kperp,xr-dx)
  do i=1,3
    do k=1,3
      ep(i,k)=ep(i,k)+chi(i,k)/(vAc*vAc*(xr-dx)*(xr-dx))
      epprime(i,k)=epprime(i,k)+chiplus(i,k)/(vAc*vAc*(xr+dx)*(xr+dx))
    enddo
  enddo
enddo


! derivative of ep at omega = omega_{kr}:
do i=1,3
 do k=1,3
  epprime(i,k)=(epprime(i,k)-ep(i,k))/(2.d0*dx)
 enddo
enddo

ep=0.d0
chi_a=0.d0

! evaluate ep at omega = omega_{kr}
do j=1,numspec
  call calc_chi(chi,j,kz,kperp,xr)
  do i=1,3
    do k=1,3
      !anti-hermitian part of the susceptibility:
  	  chi_a(j,i,k)=(chi(i,k)- conjg(chi(k,i)))/(2.d0*uniti*vAc*vAc*xr*xr)
      ep(i,k)=ep(i,k)+chi(i,k)/(vAc*vAc*xr*xr)
    enddo
  enddo
enddo


do i=1,3
  ep(i,i)=ep(i,i)+1.d0
enddo

call calc_polarization(Avec,Ek,Bk,x,kz,kperp)


! calculate Hermitian part only:
do i=1,3
	do k=1,3
		eph(i,k)=(ep(i,k)+conjg(ep(k,i)))/2.d0
		epprimeh(i,k)=(epprime(i,k)+conjg(epprime(k,i)))/2.d0
	enddo
enddo



! Now we have the dielectric tensor and the polarization. Let's calculate the wave energy
Ekmult=0.d0
Ekmultprime=0.d0
do i=1,3
  	do k=1,3
    	Ekmult(i)=Ekmult(i)+eph(i,k)*Ek(k)
    	Ekmultprime(i)=Ekmultprime(i)+epprimeh(i,k)*Ek(k)
	enddo
enddo


energy=0.d0
do i=1,3
	energy=energy+1.d0*real(conjg(Ek(i))*Ekmult(i))
	energy=energy+1.d0*real(xr*conjg(Ek(i))*Ekmultprime(i))
	energy=energy+1.d0*real(conjg(Bk(i))*Bk(i))
enddo

energy=energy/(16.d0*M_PI)


! Determine the contribution to gamma from all species (see Eq. (4) in Quataert (1998))
gamma_contribution=0.d0
heating=0.d0

do j=1,numspec
	! chi_a and Ek are already correct from above.
	do k=1,3
		do i=1,3
			heating(j) = heating(j)+ &
              1.d0*real(conjg(Ek(k))*chi_a(j,k,i)*Ek(i))
		enddo
	enddo
  heating(j)=heating(j)/(4.d0*energy)
  gamma_contribution(j)=-heating(j)*real(xr)/(4.d0*M_PI)
enddo



end subroutine
