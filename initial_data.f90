
! ===========================================================================
! initial_data.f90
! ===========================================================================
! Here are initialize all the functions defined on the grid.

!     Copyright (C) 2012  Edison Montoya, eamonto@gmail.com

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Up to date: 26 Oct 2012					


  subroutine initial_data(omega,pi,sigma,x,dx,grid_points)

    use param
    use mylibrary

    implicit none

    type(dynamical_func) omega,pi,sigma
    type(extra_func)     x,dx   

    integer :: grid_points

    real(double) :: amp    = 1.0D0
    real(double) :: x0     = 0.78539816339D0 !pi/4
    real(double) :: sigma0 = 1.0/100.0D0
    integer :: i

    dx%f = dx_aux

    do i=0,grid_points
       x%f(i) = i*dx%f(i)
    enddo

    omega%f = amp*exp(-(x%f-x0)**2/sigma0**2) + 1.0D0
    
    pi%f  = 0.0D0
    
    sigma%f = -2.0D0*(x%f-x0)*amp*exp(-(x%f-x0)**2/sigma0**2)/sigma0**2

  end subroutine initial_data
