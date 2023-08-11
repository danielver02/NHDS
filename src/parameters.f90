! This file is part of NHDS
! Copyright (C) 2020 Daniel Verscharen (d.verscharen@ucl.ac.uk)
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

subroutine set_parameters(filename)
use input_params
implicit none
character*300, intent(in) :: filename


nameList /parameters/ &
    numspec, numiter, det_D_threshold, nmax, Bessel_zero, initial_guess, scan_type, krange,&
    ksteps, alpha, beta, charge, mass, density, vdrift, theta_range, theta_steps,&
    vAc, output_mom, output_EB, kth_file, kth_filename, ampl_mode, ampl, output_warning, &
    mmax, output_df, species_df, Bessel_zero_deltaf, vxsteps, vysteps, vzsteps, vxrange, &
    vyrange, vzrange, timesteps, periods, num_periods, damping, const_r


    ! Define the standard parameters:
    numspec=2
    numiter=1000
    det_D_threshold=1.d-16
    nmax=500
    Bessel_zero=1.d-50
    initial_guess=(0.d0,0.d0)
    scan_type=1
    krange=(/0.1d0,1.d0/)
    ksteps=150
    theta_range=(/0.01d0,0.01d0/)
    theta_steps=1
    alpha=(/1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
    beta=(/1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
    charge=(/1.d0,-1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
    mass=(/1.d0,5.446623d-4,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
    density=(/1.d0,1.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
    vdrift=(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/)
    vAc=1.d-4
    output_mom=.TRUE.
    output_EB=.TRUE.
    kth_file=.FALSE.
    kth_filename=''
    ampl_mode=1
    ampl=1.d0
    output_warning=.FALSE.
    output_df=.FALSE.
    species_df=1
    mmax=1000
    Bessel_zero_deltaf=1.d-50
    vxsteps=1
    vysteps=1
    vzsteps=1
    vxrange=(/-1.d0,1.d0/)
    vyrange=(/-1.d0,1.d0/)
    vzrange=(/-1.d0,1.d0/)
    timesteps=40
    periods=.TRUE.
    num_periods=8
    damping=.FALSE.
    const_r=.TRUE.






    open (unit=5,file=trim(filename),status='old',action='read')
    read (unit=5,nml=parameters)
    close(5)

end subroutine
