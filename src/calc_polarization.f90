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