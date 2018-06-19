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

subroutine disptensor(D,x,kz,kperp)
use input_params
implicit none
double complex :: ep(3,3),D(3,3),chi(3,3)
double precision :: kz,kperp
double complex :: x
integer :: i,k,j


do i=1,3
 do k=1,3
  ep(i,k)=0.d0
 enddo
enddo



do j=1,numspec
  call calc_chi(chi,j,kz,kperp,x)
  do i=1,3
    do k=1,3
      ep(i,k)=ep(i,k)+chi(i,k)
    enddo
  enddo
enddo


!vAc is the Alfven speed in units of speed of light
do i=1,3
  ep(i,i)=ep(i,i)+vAc*vAc
enddo

do i=1,3
  do k=1,3
    D(i,k)=ep(i,k)
  enddo
enddo

D(1,1)=D(1,1)-kz*kz/(x*x)
D(1,3)=D(1,3)+kz*kperp/(x*x)
D(2,2)=D(2,2)-(kperp*kperp+kz*kz)/(x*x)
D(3,1)=D(3,1)+kz*kperp/(x*x)
D(3,3)=D(3,3)-Kperp*kperp/(x*x)

!now the tensor has been determined, the determinant has to vanish!


end subroutine
