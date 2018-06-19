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

subroutine calc_chi(chi,j,kz,kperp,x)
use input_params
implicit none
double complex :: chi(3,3),Y(3,3),Ynew(3,3),x
double precision :: kz,kperp,z,besselI
integer :: j,n,i,k,nmaxrun
logical :: Bessel_run

z=0.5d0*(kperp*vtherm(j)/Omega(j))*(kperp*vtherm(j)/Omega(j))*alpha(j)

do i=1,3
 do k=1,3
  Y(i,k)=0.d0
 enddo
enddo


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

do n=-nmaxrun,nmaxrun
   call calc_ypsilon(Ynew,j,n,kz,kperp,x)
   do i=1,3
    do k=1,3
      Y(i,k)=Y(i,k)+exp(-z)*Ynew(i,k)
    enddo
   enddo
enddo



chi(1,1)=Y(1,1)/(x*ell(j)*ell(j))
chi(1,2)=Y(1,2)/(x*ell(j)*ell(j))
chi(1,3)=Y(1,3)/(x*ell(j)*ell(j))
chi(2,1)=Y(2,1)/(x*ell(j)*ell(j))
chi(2,2)=Y(2,2)/(x*ell(j)*ell(j))
chi(2,3)=Y(2,3)/(x*ell(j)*ell(j))
chi(3,1)=Y(3,1)/(x*ell(j)*ell(j))
chi(3,2)=Y(3,2)/(x*ell(j)*ell(j))
chi(3,3)=2.d0*vdrift(j)/(ell(j)*ell(j)*x*kz*vtherm(j)*vtherm(j)*alpha(j))+Y(3,3)/(x*ell(j)*ell(j))



end subroutine







subroutine calc_ypsilon(Y,j,n,kz,kperp,x)
use input_params
implicit none
double complex :: zeta,x,Y(3,3),An,Bn,resfac,dispfunct
double precision :: kz,kperp,BInz,z
double precision :: besselI
double precision :: dBInzdz
integer :: n,j
logical :: kpos

kpos=.TRUE.
if (kz.LT.0.d0) kpos=.FALSE.



zeta=(x-kz*vdrift(j)-1.d0*n*Omega(j))/(kz*vtherm(j))
resfac=x-kz*vdrift(j)-1.d0*n*Omega(j)
z=0.5d0*(kperp*vtherm(j)/Omega(j))*(kperp*vtherm(j)/Omega(j))*alpha(j)

An=(1.d0/x)*(alpha(j)-1.d0)
An=An+(1.d0/(kz*x*vtherm(j))) *( alpha(j)*resfac + 1.d0*n*Omega(j))*dispfunct(zeta,kpos)

Bn=(alpha(j)*(x-1.d0*n*Omega(j))-(kz*vdrift(j)-1.d0*n*Omega(j)))/(kz*x)
Bn=Bn+((x-1.d0*n*Omega(j))*(alpha(j)*resfac+1.d0*n*Omega(j))/(kz*kz*vtherm(j)*x) )*dispfunct(zeta,kpos)

if (n.GE.0) then
	BInz=1.d0*besselI(n,z)
	dBInzdz=besselI(n+1,z)+1.d0*n*BInz/z
else
	BInz=1.d0*besselI(-n,z)
	dBInzdz=besselI(-n-1,z)+1.d0*n*BInz/z
endif


! The tensor in Stix's (10-57)
Y(1,1)=1.d0*(n*n)*BInz*An/z
Y(1,2)=-uniti*n*(BInz-dBInzdz)*An
Y(1,3)=kperp*n*BInz*Bn/(Omega(j)*z)
Y(2,1)=uniti*n*(BInz-dBInzdz)*An
Y(2,2)=(1.d0*(n*n)*BInz/z+2.d0*z*BInz-2.d0*z*dBInzdz)*An
Y(2,3)=uniti*kperp*(BInz-dBInzdz)*Bn/Omega(j)
Y(3,1)=kperp*BInz*n*Bn/(Omega(j)*z)
Y(3,2)=-uniti*kperp*(BInz-dBInzdz)*Bn/Omega(j)
Y(3,3)=2.d0*(x-1.d0*n*Omega(j))*BInz*Bn/(kz*vtherm(j)*vtherm(j)*alpha(j))


end subroutine





function besselI(n,x)
implicit none
double precision :: besselI,x,BESSI
integer :: n

if (n.LT.0) then
	besselI=BESSI(-n,x)
else
	besselI=BESSI(n,x)
endif

end function



!************************************************************************
!*                                                                      *
!*    Program to calculate the first kind modified Bessel function of   *
!*    integer order N, for any REAL X, using the function BESSI(N,X).   *
!*                                                                      *
!* -------------------------------------------------------------------- *
!* -------------------------------------------------------------------- *
!*   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
!*                                                                      *
!*                               F90 Release 1.1 By J-P Moreau, Paris.  *
!*                                                                      *
!*   Version 1.1: corected value of P4 in BESSIO (P4=1.2067492 and not  *
!*                1.2067429) Aug. 2011.                                 *
!************************************************************************


! ----------------------------------------------------------------------
      FUNCTION BESSI(N,X)
!
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows.
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!
      PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
      REAL *8 X,BESSI,BESSI0,BESSI1,TOX,BIM,BI,BIP
      IF (N.EQ.0) THEN
      BESSI = BESSI0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSI = BESSI1(X)
      RETURN
      ENDIF
      IF(X.EQ.0.D0) THEN
      BESSI=0.D0
      RETURN
      ENDIF
      TOX = 2.D0/X
      BIP = 0.D0
      BI  = 1.D0
      BESSI = 0.D0
      M = 2*((N+INT(SQRT(FLOAT(IACC*N)))))
      DO 12 J = M,1,-1
      BIM = BIP+DFLOAT(J)*TOX*BI
      BIP = BI
      BI  = BIM
      IF (ABS(BI).GT.BIGNO) THEN
      BI  = BI*BIGNI
      BIP = BIP*BIGNI
      BESSI = BESSI*BIGNI
      ENDIF
      IF (J.EQ.N) BESSI = BIP
   12 CONTINUE
      BESSI = BESSI*BESSI0(X)/BI
      RETURN
      END
! ----------------------------------------------------------------------
! Auxiliary Bessel functions for N=0, N=1
      FUNCTION BESSI0(X)
      REAL *8 X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,  &
      Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067492D0,  &
      0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
      0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
      0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      ENDIF
      RETURN
      END
! ----------------------------------------------------------------------
      FUNCTION BESSI1(X)
      REAL *8 X,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,  &
      Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
      0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
      -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1, &
      -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
      ENDIF
      RETURN
      END
! ----------------------------------------------------------------------
