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
