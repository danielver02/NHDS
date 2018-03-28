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

subroutine calc_xi(xi,j,pol,polz,x,kz,kperp)
use input_params
implicit none
double precision :: kz,kperp
double complex :: x,pol,polz,Avec(3),xi,chi(3,3)
integer :: j

! This subroutine calculates the quantity xi in dn/n_0 = xi * dB_(x,y,z)/B_0


if (ampl_mode.EQ.1) then
	! Avec fulfills the equation  c * dE = Avec * vA * dB_x
	Avec(1)=uniti*x/(kz*pol)
	Avec(2)=-x/kz
	Avec(3)=polz*uniti*x/(kz*pol)	
else if (ampl_mode.EQ.2) then
	! Avec fulfills the equation  c * dE = Avec * vA * dB_y
	Avec(1)=x/(kz-polz*kperp)
	Avec(2)=uniti*pol*x/(kz-polz*kperp)
	Avec(3)=x/((kz/polz)-kperp)
else if (ampl_mode.EQ.3) then
	! Avec fulfills the equation  c * dE = Avec * vA * dB_z
	Avec(1)=-uniti*x/(kperp*pol)
	Avec(2)=x/kperp
	Avec(3)=-uniti*x*polz/(kperp*pol)
endif


call calc_chi(chi,j,kz,kperp,x)

! continuity equation (as in Stix 10-75)
xi=kperp*(chi(1,1)*Avec(1)+chi(1,2)*Avec(2)+chi(1,3)*Avec(3))
xi=xi+kz*(chi(3,1)*Avec(1)+chi(3,2)*Avec(2)+chi(3,3)*Avec(3))
xi=-uniti*xi*ampl/(density(j)*charge(j))

! Normalization: chi is different by a factor vAc*vAc from Stix' definition:
! our chi divided by (v_A/c)^2 is equal to Stix' chi.

end subroutine



subroutine calc_fluctRj(deltaRj,dUperpx,dUperpy,dUpar,j,pol,polz,x,kz,kperp)
use input_params
implicit none
double precision :: kz,kperp,z,besselI,BInz,dBInzdz,d2BInzdz2
double complex :: x,pol,polz,Avec(3),dispfunct,zeta,deltaRj,chi(3,3)
double complex :: Z0,Z1,Z2,Z3,dpperp,dppar,Exterm,Eyterm,Ezterm,dUpar
double complex :: dUperpx,dUperpy,Ubar
logical :: kpos,Bessel_run
integer :: j,n,nmaxrun

! This subroutine calculates the pressure fluctuations that are
! needed to calculate the quantity R_j (the fluctuations in the temperature anisotropy of species j)
! define ampl = dBz/B0

z=0.5d0*(kperp*vtherm(j)/Omega(j))*(kperp*vtherm(j)/Omega(j))*alpha(j)

kpos=.TRUE.
if (kz.LT.0.d0) kpos=.FALSE.


if (ampl_mode.EQ.1) then
	! Avec fulfills the equation  c * dE = Avec * vA * dB_x
	Avec(1)=uniti*x/(kz*pol)
	Avec(2)=-x/kz
	Avec(3)=polz*uniti*x/(kz*pol)	
else if (ampl_mode.EQ.2) then
	! Avec fulfills the equation  c * dE = Avec * vA * dB_y
	Avec(1)=x/(kz-polz*kperp)
	Avec(2)=uniti*pol*x/(kz-polz*kperp)
	Avec(3)=x/((kz/polz)-kperp)
else if (ampl_mode.EQ.3) then
	! Avec fulfills the equation  c * dE = Avec * vA * dB_z
	Avec(1)=-uniti*x/(kperp*pol)
	Avec(2)=x/kperp
	Avec(3)=-uniti*x*polz/(kperp*pol)
endif

! This is for the calculation of dUpar and dUperp from the definition of the current density
! and the susceptibility. These expressions are equivalent to taking the first moment
! of delta f and lead to the same result.
call calc_chi(chi,j,kz,kperp,x)
dUpar=chi(3,1)*Avec(1)+chi(3,2)*Avec(2)+chi(3,3)*Avec(3)
dUpar=-uniti*x*dUpar*ampl/(density(j)*charge(j))

dUperpx=chi(1,1)*Avec(1)+chi(1,2)*Avec(2)+chi(1,3)*Avec(3)
dUperpx=-uniti*x*dUperpx*ampl/(density(j)*charge(j))

dUperpy=chi(2,1)*Avec(1)+chi(2,2)*Avec(2)+chi(2,3)*Avec(3)
dUperpy=-uniti*x*dUperpy*ampl/(density(j)*charge(j))



! This is without dUpar:
Ubar=vdrift(j)!+dUpar

! Normalization: chi is different by a factor vAc*vAc from Stix' definition:
! our chi divided by (v_A/c)^2 is equal to Stix' chi.

! dUperpx and dUperpy need to be small so that the following calculation is accurate.

!determine maximum n for Besselfunction:
nmaxrun=nmax
n=0
Bessel_run=.TRUE.
do while (Bessel_run)
  if ((n.GE.nmax).OR.(besselI(n,z).LT.Bessel_zero)) then
     nmaxrun=n
     Bessel_run=.FALSE.
  endif
  n=n+1
enddo

if ((n.GE.nmax).AND.(output_warning)) then
 write (*,*) "WARNING: Insufficient order of Bessel function (nmax)."
 stop
endif


! dpperp is deltapperp^2 / ampl in units of the magnetic background pressure
! dppar is deltappar^2 / ampl in units of the magnetic background pressure.

dpperp=0.d0
dppar=0.d0


! Run over all n:
do n=-nmaxrun,nmaxrun
   
	if (n.GE.0) then
		BInz=1.d0*besselI(n,z)
		dBInzdz=besselI(n+1,z)+1.d0*n*BInz/z 
	else
		BInz=1.d0*besselI(-n,z)
		dBInzdz=besselI(-n-1,z)+1.d0*n*BInz/z
	endif

	if (n.GE.2) d2BInzdz2=0.25d0*(besselI(n-2,z)+2.d0*besselI(n,z)+besselI(n+2,z))
	if (n.LE.-2) d2BInzdz2=0.25d0*(besselI(2-n,z)+2.d0*besselI(-n,z)+besselI(-2-n,z))
	if (n.EQ.1) d2BInzdz2=0.25d0*(3.d0*besselI(1,z)+besselI(3,z))
	if (n.EQ.-1) d2BInzdz2=0.25d0*(besselI(3,z)+3.d0*besselI(1,z))
	if (n.EQ.0) d2BInzdz2=0.5d0*(besselI(2,z)+besselI(0,z))
	
	zeta=(x-kz*vdrift(j)-1.d0*n*Omega(j))/(kz*vtherm(j))

	! determine Stix' (normalized) dispersion functions:
	! (normalized as follows: Z0, Z1/vA, Z2*vA^2, Z3*vA^3)
	Z0=dispfunct(zeta,kpos)
	Z1=vtherm(j)*(1.d0+(zeta+vdrift(j)/vtherm(j))*Z0)
	Z2=(zeta+2.d0*vdrift(j)/vtherm(j)+Z0*(zeta+vdrift(j)/vtherm(j))**(2.d0))
	Z2=Z2*vtherm(j)*vtherm(j)
	Z3=0.5d0+zeta*zeta+3.d0*zeta*vdrift(j)/vtherm(j)+3.d0*(vdrift(j)/vtherm(j))**(2.d0)
	Z3=Z3+Z0*(zeta+vdrift(j)/vtherm(j))**(3.d0)
	Z3=Z3*vtherm(j)*vtherm(j)*vtherm(j)


	! This is for the calculation of dpperp:
	
	Exterm=Avec(1)*(BInz+z*(dBinzdz-BInz))*(1.d0*n)*Omega(j)/kperp
	Eyterm=-uniti*Avec(2)*(BInz-dBInzdz-0.5d0*z*(BInz+d2BInzdz2)+z*dBInzdz)
	Eyterm=Eyterm*kperp*vtherm(j)*vtherm(j)*alpha(j)/Omega(j)
	Ezterm=(1.d0-1.d0*n*Omega(j)/x)*alpha(j)*(Z1-vdrift(j)*Z0)+1.d0*n*Omega(j)*Z1/x
	Ezterm=Ezterm*Avec(3)*(BInz+z*(dBInzdz-BInz))
	
	dpperp=dpperp+(Exterm+Eyterm)*(Z0-kz*Z1/x+alpha(j)*kz*(Z1-vdrift(j)*Z0)/x)+Ezterm

	
	! This is for the calculation of dppar:
	
	Exterm=alpha(j)*kz*(Z3-vdrift(j)*Z2-2.d0*Ubar*Z2+2.d0*vdrift(j)*Ubar*Z1+Ubar*Ubar*Z1-vdrift(j)*Ubar*Ubar*Z0)/x
	Exterm=Exterm+Z2-2.d0*Ubar*Z1+Ubar*Ubar*Z0
	Exterm=Exterm-kz*(Z3-2.d0*Ubar*Z2+Ubar*Ubar*Z1)/x
	Exterm=Exterm*Avec(1)*BInz*(1.d0*n)*Omega(j)/(kperp*vtherm(j)*vtherm(j)*alpha(j))
	Eyterm=alpha(j)*kz*(Z3-vdrift(j)*Z2-2.d0*Ubar*Z2+2.d0*vdrift(j)*Ubar*Z1+Ubar*Ubar*Z1-vdrift(j)*Ubar*Ubar*Z0)/x
	Eyterm=Eyterm+Z2-2.d0*Ubar*Z1+Ubar*Ubar*Z0
	Eyterm=Eyterm-kz*(Z3-2.d0*Ubar*Z2+Ubar*Ubar*Z1)/x
	Eyterm=-uniti*Eyterm*Avec(2)*(BInz-dBInzdz)*kperp/(2.d0*Omega(j))
	Ezterm=(1.d0-(1.d0*n)*Omega(j)/x)*(Z3-vdrift(j)*Z2-2.d0*Ubar*Z2+2.d0*vdrift(j)*Ubar*Z1+Ubar*Ubar*Z1-vdrift(j)*Ubar*Ubar*Z0)
	Ezterm=Ezterm+(1.d0*n)*Omega(j)*(Z3-2.d0*Ubar*Z2+Ubar*Ubar*Z1)/(alpha(j)*x)
	Ezterm=Ezterm*Avec(3)*BInz/(vtherm(j)*vtherm(j))



	dppar=dppar+Exterm+Eyterm+Ezterm

enddo


dpperp=-2.d0*uniti*mass(j)*density(j)*Omega(j)*dpperp*exp(-z)/(kz*vtherm(j))
dppar=-4.d0*uniti*mass(j)*density(j)*Omega(j)*dppar*exp(-z)/(kz*vtherm(j))


deltaRj = (vtherm(j)*vtherm(j)*alpha(j)+ampl*dpperp/(density(j)*mass(j)))
deltaRj = deltaRj /(vtherm(j)*vtherm(j)+ampl*dppar/(density(j)*mass(j)))

end subroutine

