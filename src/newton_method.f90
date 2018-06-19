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
