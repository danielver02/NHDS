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

subroutine calc_polarization(pol,polz,x,kz,kperp)
use input_params
implicit none

double complex :: pol, polz, D(3,3),x
double precision :: kz,kperp


call disptensor(D,x,kz,kperp)

pol=(D(1,1)*D(2,3)-D(1,3)*D(2,1))/(D(1,3)*D(2,2)-D(1,2)*D(2,3))
pol=pol/uniti


! pol is the quantity -iEky/Ekx = (Ekr-Ekl)/(Ekl+Ekr)

polz=(D(1,1)*D(2,2)-D(1,2)*D(2,1))/(D(1,2)*D(2,3)-D(1,3)*D(2,2))

! polz is the quantity Ekz/Ekx

end subroutine
