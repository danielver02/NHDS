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