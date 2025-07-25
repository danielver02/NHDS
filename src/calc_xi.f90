! This file is part of NHDS
! Copyright (C) 2025 Daniel Verscharen (d.verscharen@ucl.ac.uk)
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



subroutine calc_xi(xi,dpperp,dppar,dUperpx,dUperpy,dUpar,j,Avec,x,kz,kperp)
use input_params
implicit none
double precision :: kz,kperp,z,besselI,BInz,dBInzdz,d2BInzdz2
double complex :: x,Avec(3),dispfunct,zeta,chi(3,3),xi
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

 
! Normalization: chi is different by a factor vAc*vAc*x*x from Stix' definition:
! our chi divided by (v_A/c)^2 and x^2 is equal to Stix' chi.
call calc_chi(chi,j,kz,kperp,x)

! continuity equation (as in Stix 10-75)
! The following calculates the quantity xi in dn/n_0 = xi * dB_(x,y,z)/B_0
xi=kperp*(chi(1,1)*Avec(1)+chi(1,2)*Avec(2)+chi(1,3)*Avec(3))
xi=xi+kz*(chi(3,1)*Avec(1)+chi(3,2)*Avec(2)+chi(3,3)*Avec(3))
xi=-uniti*xi*ampl/(x*x*density(j)*charge(j))


! The following calculates dUpar and dUperp from the definition of the current density
! and the susceptibility. These expressions are equivalent to taking the first moment
! of delta f and lead to the same result.

dUpar=chi(3,1)*Avec(1)+chi(3,2)*Avec(2)+chi(3,3)*Avec(3)
dUpar=-uniti*dUpar*ampl/(x*density(j)*charge(j))
dUpar=dUpar-xi*vdrift(j) ! This final term accounts for the parallel flow of the species j.

dUperpx=chi(1,1)*Avec(1)+chi(1,2)*Avec(2)+chi(1,3)*Avec(3)
dUperpx=-uniti*dUperpx*ampl/(x*density(j)*charge(j))

dUperpy=chi(2,1)*Avec(1)+chi(2,2)*Avec(2)+chi(2,3)*Avec(3)
dUperpy=-uniti*dUperpy*ampl/(x*density(j)*charge(j))



! Check whether we run the cold plasma case:
if (vtherm(j).EQ.0.d0) then

  dpperp=0.d0
  dppar=0.d0

else


! This is without dUpar:
Ubar=vdrift(j)!+dUpar


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


! dpperp is deltapperp in units of the magnetic background pressure
! dppar is deltappar in units of the magnetic background pressure.

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

! Remember that the Bessel functions I_n already include a factor exp(-z), so it does not appear here anymore:
dpperp=-2.d0*uniti*ampl*mass(j)*density(j)*Omega(j)*dpperp/(kz*vtherm(j))
dppar=-4.d0*uniti*ampl*mass(j)*density(j)*Omega(j)*dppar/(kz*vtherm(j))

endif


end subroutine
