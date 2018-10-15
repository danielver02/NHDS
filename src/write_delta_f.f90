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

subroutine write_delta_f(j,kk,x,pol,polz)
use input_params
use HDF5
implicit none
double complex :: deltaf,comp1,comp2,comp3,comp4,x,pol,polz,Avec(3),UStix,VStix,a
double precision :: phi,vperp,vpar,dfdvpar,dfdvperp,kk,kz,kperp,z,bessel,fnull,vx,vy,vz
double precision :: time, theta
real, allocatable :: VXarray(:,:,:),VYarray(:,:,:),VZarray(:,:,:),farray(:,:,:)
real :: minf,maxf,maxdfim,mindfim
integer :: i,k,l,j,m,timerun,mmaxrun,error
logical :: Bessel_run
integer(HID_T) :: file_id		! File identifier
integer(HID_T) :: dset1_id,dset2_id,dset3_id,dset4_id       ! Dataset identifier
integer(HID_T) :: dspace1_id,dspace2_id,dspace3_id,dspace4_id      ! Dataspace identifier
integer(HSIZE_T) :: dims(3) ! Dataset dimensions
integer ::   rank=3 ! Dataset rank
character(LEN=19) :: filename ! File name for the HDF5 file
character(LEN=3) :: timestring
character(LEN=2), PARAMETER :: dset1name = "VX"     ! Dataset name
character(LEN=2), PARAMETER :: dset2name = "VY"     ! Dataset name
character(LEN=2), PARAMETER :: dset3name = "VZ"     ! Dataset name
character(LEN=6), PARAMETER :: dset4name = "DELTAF"     ! Dataset name

dims=(/ vxsteps+1,vysteps+1,vzsteps+1 /)

allocate(VXarray(vxsteps+1,vysteps+1,vzsteps+1))
allocate(VYarray(vxsteps+1,vysteps+1,vzsteps+1))
allocate(VZarray(vxsteps+1,vysteps+1,vzsteps+1))
allocate(farray(vxsteps+1,vysteps+1,vzsteps+1))


kz=kk*cos(theta*M_PI/180.d0)
kperp=kk*sin(theta*M_PI/180.d0)



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

write (*,*) "# Writing distribution function to files output_deltafXXX.h5"
call h5open_f(error)

! Start Time series:
maxf=-999.
minf=999.
maxdfim=-999.
mindfim=999.




do timerun=0,num_periods*timesteps


if (periods) then
	time=2.d0*M_PI*(1.d0*timerun)/(real(x)*1.d0*timesteps)
else
	time=2.d0*M_PI*(1.d0*timerun)/(real(1.d0)*1.d0*timesteps)
endif



filename="output_deltaf000.h5"
if (timerun.LT.10) then
	 write ( timestring, '(i1)' ) timerun
	 filename(16:17)=timestring
	 filename(17:19)='.h5'
else if ((timerun.GE.10).AND.(timerun.LT.100)) then
	 write ( timestring, '(i2)' ) timerun
	 filename(15:16)=timestring
else if ((timerun.GE.100).AND.(timerun.LT.1000)) then
	 write ( timestring, '(i3)' ) timerun
	filename(14:16)=timestring
else if (timerun.GE.1000) then
	write (*,*) "ERROR: More than 999 timesteps."
endif

! Initialize HDF5 interface and create file:
call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

! Create the dataspace and dataset with default properties:
call h5screate_simple_f(rank, dims, dspace1_id, error) ! for VX
call h5dcreate_f(file_id, dset1name, H5T_NATIVE_REAL, dspace1_id,dset1_id, error) ! for VX
call h5screate_simple_f(rank, dims, dspace2_id, error) ! for VY
call h5dcreate_f(file_id, dset2name, H5T_NATIVE_REAL, dspace2_id,dset2_id, error) ! for VY
call h5screate_simple_f(rank, dims, dspace3_id, error) ! for VZ
call h5dcreate_f(file_id, dset3name, H5T_NATIVE_REAL, dspace3_id,dset3_id, error) ! for VZ
call h5screate_simple_f(rank, dims, dspace4_id, error)
call h5dcreate_f(file_id, dset4name, H5T_NATIVE_REAL, dspace4_id,dset4_id, error) ! for distr. function


! Here should the loop over all v start:

do i=0,vxsteps
do k=0,vysteps
do l=0,vzsteps

vx=vxrange(1)+(1.d0*i)*(vxrange(2)-vxrange(1))/(1.d0*vxsteps)
vy=vyrange(1)+(1.d0*k)*(vyrange(2)-vyrange(1))/(1.d0*vysteps)
vz=vzrange(1)+(1.d0*l)*(vzrange(2)-vzrange(1))/(1.d0*vzsteps)

vpar=vz
vperp=sqrt(vx*vx+vy*vy)
phi=acos(vx/vperp)
if (vy.LT.0.d0) phi=2.d0*M_PI-phi
if (vperp.EQ.0.d0) phi=0.d0

dfdvpar=-2.d0*(vpar-vdrift(j))*exp(-vperp*vperp/(alpha(j)*vtherm(j)*vtherm(j)))
dfdvpar=dfdvpar*exp(-(vpar-vdrift(j))*(vpar-vdrift(j))/(vtherm(j)*vtherm(j)))
dfdvpar=dfdvpar/(M_PI*sqrt(M_PI)*alpha(j)*(vtherm(j)**5.d0))

dfdvperp=-2.d0*vperp*exp(-vperp*vperp/(alpha(j)*vtherm(j)*vtherm(j)))
dfdvperp=dfdvperp*exp(-(vpar-vdrift(j))*(vpar-vdrift(j))/(vtherm(j)*vtherm(j)))
dfdvperp=dfdvperp/(M_PI*sqrt(M_PI)*alpha(j)*alpha(j)*(vtherm(j)**5.d0))

fnull=1.d0/(alpha(j)*(vtherm(j)**3.d0)*(M_PI)**(3.d0/2.d0))
fnull=fnull*exp(-(vpar-vdrift(j))*(vpar-vdrift(j))/(vtherm(j)*vtherm(j)))
fnull=fnull*exp(-vperp*vperp/(alpha(j)*vtherm(j)*vtherm(j)))

z=kperp*vperp/Omega(j)


!determine maximum m for Besselfunction:
mmaxrun=mmax
m=0
Bessel_run=.TRUE.
do while (Bessel_run)
  if ((m.GE.mmax).OR.(bessel(m,z).LT.Bessel_zero_deltaf)) then
     mmaxrun=m
     Bessel_run=.FALSE.
  endif
  m=m+1
enddo


if ((m.GE.mmax).AND.(output_warning)) then
 write (*,*) "WARNING: Insufficient order of Bessel function (mmax) in write_delta_f.f90."
 stop
endif

!Summation over all Bessel functions:
deltaf=0.d0
do m=-mmaxrun,mmaxrun

	a=(1.d0*m)*Omega(j)-x+kz*vpar


	UStix=dfdvperp+kz*(vperp*dfdvpar-vpar*dfdvperp)/x
	VStix=kperp*(vperp*dfdvpar+vpar*dfdvperp)/x

	comp1=Avec(1)*Ustix*(a*uniti*cos(phi)-Omega(j)*sin(phi))/(Omega(j)*Omega(j)-a*a)
	comp2=Avec(2)*UStix*(a*uniti*sin(phi)-Omega(j)*cos(phi))/(Omega(j)*Omega(j)-a*a)
	comp3=-uniti*Avec(3)*dfdvpar/a
	comp4=-Avec(3)*VStix*(a*uniti*cos(phi)-Omega(j)*sin(phi))/(Omega(j)*Omega(j)-a*a)

	deltaf = deltaf - exp(uniti*z*sin(phi))*exp(-uniti*(1.d0*m)*phi)*bessel(m,z)*ampl*(comp1+comp2+comp3+comp4)

enddo ! end m loop


if (const_r) then
	if (damping) then
		deltaf=deltaf*exp(-uniti*time*x)
	else
		deltaf=deltaf*exp(-uniti*time*real(x))
	endif
else
	deltaf=deltaf*exp(uniti*2.d0*M_PI*(1.d0*timerun)/(1.d0*timesteps))
endif




maxf=max((maxf),real(fnull+deltaf))
minf=min((minf),real(fnull+deltaf))
maxdfim=max((maxdfim),aimag(deltaf))
mindfim=min((mindfim),aimag(deltaf))

	VXarray(i+1,k+1,l+1)=real(vx)
	VYarray(i+1,k+1,l+1)=real(vy)
	VZarray(i+1,k+1,l+1)=real(vz)
	farray(i+1,k+1,l+1)=real(fnull+deltaf)


enddo ! end vz loop
enddo ! end vy loop
enddo ! end vx loop

call h5dwrite_f(dset1_id, H5T_NATIVE_REAL, VXarray, dims, error)
call h5dwrite_f(dset2_id, H5T_NATIVE_REAL, VYarray, dims, error)
call h5dwrite_f(dset3_id, H5T_NATIVE_REAL, VZarray, dims, error)
call h5dwrite_f(dset4_id, H5T_NATIVE_REAL, farray, dims, error)


call h5dclose_f(dset1_id, error)
call h5dclose_f(dset2_id, error)
call h5dclose_f(dset3_id, error)
call h5dclose_f(dset4_id, error)
call h5sclose_f(dspace1_id, error)
call h5sclose_f(dspace2_id, error)
call h5sclose_f(dspace3_id, error)
call h5sclose_f(dspace4_id, error)



! Close HDF5 interface:
call h5fclose_f(file_id, error)
enddo
call h5close_f(error)

call generate_xdmf

write (*,*) "# Maximum value of f:",maxf
write (*,*) "# Minimum value of f:",minf
write (*,*) "# Maximum value of Im(delta f):",maxdfim
write (*,*) "# Minimum value of Im(delta f):",mindfim

end subroutine


subroutine generate_xdmf
use input_params
implicit none
integer :: i
character(LEN=14) :: form
character(LEN=1) :: vxstring,vystring,vzstring

write (*,*) "# Generating XDMF file visualize_output_deltaf.xmf"
open(unit=100,file='visualize_output_deltaf.xmf',status='replace',action='write')



!define the filenames for the output
  write ( vxstring, '(i1)' ) int(log10(1.*vxsteps))+2
  write ( vystring, '(i1)' ) int(log10(1.*vysteps))+2
  write ( vzstring, '(i1)' ) int(log10(1.*vzsteps))+1

  form(1:4)="(A,I"
  form(5:6)=vzstring
  form(6:7)=",I"
  form(8:9)=vystring
  form(9:10)=",I"
  form(11:12)=vxstring
  form(12:14)=",A)"


write (100,'(A)') '<?xml version="1.0" ?>'
write (100,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
write (100,'(A)') '<Xdmf Version="2.0">'
write (100,'(A)') '<Domain>'
write (100,'(A)') '<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'

do i=0,num_periods*timesteps
	write (100,'(A)') '  <Grid Name="velocity_space" GridType="Uniform">'
	write (100,'(A,I3,A)') '  <Time Type="Single" Value="',i,'" />'
	write (100,form) '    <Topology TopologyType="3DSMesh" NumberOfElements="', vzsteps+1,vysteps+1,vxsteps+1,'"/>'
	write (100,'(A)') '     <Geometry GeometryType="X_Y_Z">'
	write (100,form) '       <DataItem Dimensions="', vzsteps+1,vysteps+1,vxsteps+1,'" NumberType="Float" Precision="4" Format="HDF">'
	write (100,'(A)') '        output_deltaf001.h5:/VX'
	write (100,'(A)') '       </DataItem>'
	write (100,form) '       <DataItem Dimensions="', vzsteps+1,vysteps+1,vxsteps+1,'" NumberType="Float" Precision="4" Format="HDF">'
	write (100,'(A)') '        output_deltaf001.h5:/VY'
	write (100,'(A)') '       </DataItem>'
	write (100,form) '        <DataItem Dimensions="', vzsteps+1,vysteps+1,vxsteps+1,'" NumberType="Float" Precision="4" Format="HDF">'
	write (100,'(A)') '        output_deltaf001.h5:/VZ'
	write (100,'(A)') '       </DataItem>'
	write (100,'(A)') '     </Geometry>'
	write (100,'(A)') '     <Attribute Name="distribution_function" AttributeType="Scalar" Center="Node">'
	write (100,form) '       <DataItem Dimensions="', vzsteps+1,vysteps+1,vxsteps+1,'" NumberType="Float" Precision="4" Format="HDF">'
	if (i.LT.10) then
		write (100,'(A,I1,A)') '        output_deltaf00',i,'.h5:/DELTAF'
	else if ((i.GE.10).AND.(i.LT.100)) then
		write (100,'(A,I2,A)') '        output_deltaf0',i,'.h5:/DELTAF'
	else if ((i.GE.100).AND.(i.LT.1000)) then
		write (100,'(A,I3,A)') '        output_deltaf',i,'.h5:/DELTAF'
	endif
	write (100,'(A)') '       </DataItem>'
	write (100,'(A)') '     </Attribute>'
	write (100,'(A)') '   </Grid>'

enddo


write (100,'(A)') '</Grid>'
write (100,'(A)') ' </Domain>'
write (100,'(A)') '</Xdmf>'


close(100)

end subroutine

function bessel(n,x)
implicit none
integer :: n
	double precision :: bessel,x,bessj

	if (n.LT.0) then
		bessel = (-1.d0**(1.d0*n))*bessj(-n,x)
	else
		bessel = bessj(n,x)
	endif

end function



     FUNCTION BESSJ (N,X)

!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows.
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.

      IMPLICIT NONE
      INTEGER, PARAMETER :: IACC = 40
      REAL*8, PARAMETER :: BIGNO = 1.D10, BIGNI = 1.D-10
      INTEGER M, N, J, JSUM
      REAL *8 X,BESSJ,BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM
      IF (N.EQ.0) THEN
      BESSJ = BESSJ0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSJ = BESSJ1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.) THEN
      BESSJ = 0.
      RETURN
      ENDIF
      TOX = 2./X
      IF (X.GT.FLOAT(N)) THEN
      BJM = BESSJ0(X)
      BJ  = BESSJ1(X)
      DO 11 J = 1,N-1
      BJP = J*TOX*BJ-BJM
      BJM = BJ
      BJ  = BJP
   11 CONTINUE
      BESSJ = BJ
      ELSE
      M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
      BESSJ = 0.
      JSUM = 0
      SUM = 0.
      BJP = 0.
      BJ  = 1.
      DO 12 J = M,1,-1
      BJM = J*TOX*BJ-BJP
      BJP = BJ
      BJ  = BJM
      IF (ABS(BJ).GT.BIGNO) THEN
      BJ  = BJ*BIGNI
      BJP = BJP*BIGNI
      BESSJ = BESSJ*BIGNI
      SUM = SUM*BIGNI
      ENDIF
      IF (JSUM.NE.0) SUM = SUM+BJ
      JSUM = 1-JSUM
      IF (J.EQ.N) BESSJ = BJP
   12 CONTINUE
      SUM = 2.*SUM-BJ
      BESSJ = BESSJ/SUM
      ENDIF
      RETURN
      END

      FUNCTION BESSJ0 (X)
      IMPLICIT NONE
      REAL *8 X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX

!     This subroutine calculates the First Kind Bessel Function of
!     order 0, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.

      REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
      DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
      DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
      IF(X.EQ.0.D0) GO TO 1
      AX = ABS (X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ0 = FR/FS
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-.785398164
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
      ENDIF
      RETURN
    1 BESSJ0 = 1.D0
      RETURN
      END
! ---------------------------------------------------------------------------
      FUNCTION BESSJ1 (X)
      IMPLICIT NONE
      REAL *8 X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
!     This subroutine calculates the First Kind Bessel Function of
!     order 1, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
      REAL *8 Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
      DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, &
      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
      DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

      AX = ABS(X)
      IF (AX.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ1 = X*(FR/FS)
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-2.35619491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
      ENDIF
      RETURN
      END
