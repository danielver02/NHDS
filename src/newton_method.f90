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

subroutine newton_method(x,kz,kperp,quality)
use input_params
implicit none
double  complex :: determinant,prevx,x,detprevx,detx,jump
integer :: iter
double precision :: kperp,kz,dblabs
logical :: go_for_newton
real :: quality

prevx=x-1.d-5-1.d-5*uniti
! newton iteration for the secant method
detprevx=determinant(prevx,kz,kperp)
iter = 0
go_for_newton=.TRUE.
do while((iter.LT.numiter).AND.(go_for_newton))
   iter=iter+1
   detx=determinant(x,kz,kperp)
   if ((dblabs(detx-detprevx).LT.1.d-80).OR.(dblabs(detx).LT.det_D_threshold)) then
        jump=0.d0
        go_for_newton=.FALSE.
   else
      jump=detx*(x-prevx)/(detx-detprevx) 
   endif     
   prevx=x
   x=x-jump
   detprevx=detx
enddo
if ((iter.GE.numiter).AND.(output_warning)) then
  write (*,*) "WARNING: Insufficient number of iterations (numiter)."
  stop
endif

quality=1.*dblabs(detx)

end subroutine










double complex function determinant(x,kz,kperp)
use input_params
implicit none
double complex :: x,D(3,3)
double precision :: kz,kperp

call disptensor(D,x,kz,kperp)

determinant=D(1,1)*D(2,2)*D(3,3)
determinant=determinant+D(1,2)*D(2,3)*D(3,1)+D(1,3)*D(2,1)*D(3,2)
determinant=determinant-D(3,1)*D(2,2)*D(1,3)-D(3,2)*D(2,3)*D(1,1)-D(3,3)*D(2,1)*D(1,2)


end function


double precision function dblabs(z)
implicit none
double complex :: z
double precision :: realpart,imagpart

realpart=1.d0*real(z)
imagpart=1.d0*aimag(z)

dblabs=dsqrt(realpart*realpart+imagpart*imagpart)

end function