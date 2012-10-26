
! ===========================================================================
! boundaries.f90
! ===========================================================================
! Implementation of the boundaries conditions

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


  subroutine boundaries(omega,pi,sigma,x,dx,grid_points)

    use mylibrary
    use param

    implicit none

    type(dynamical_func) omega,pi,sigma
    type(extra_func)     x,dx   

    integer :: grid_points

    real(double) :: der_pi,der_psi

    if(boundary.eq.1) then

       der_pi  = (- pi%f(2)+4.0D0* pi%f(1)-3.0D0* pi%f(0))/(2.0D0*dx%f(0))

       omega%f(0) = 1.0D0

       pi%s(0)  = 0.0D0
       sigma%s(0) = der_pi

       der_pi  = ( pi%f(grid_points-2)-4.0D0* pi%f(grid_points-1)+3.0D0 *pi%f(grid_points))/(2.0D0*dx%f(grid_points))

       omega%f(grid_points) = 1.0D0

       pi%s(grid_points)  = 0.0D0
       sigma%s(grid_points) = der_pi

    else if(boundary.eq.2) then

       der_pi  = (- pi%f(2)+4.0D0* pi%f(1)-3.0D0* pi%f(0))/(2.0D0*dx%f(0))
       der_psi = (-sigma%f(2)+4.0D0*sigma%f(1)-3.0D0*sigma%f(0))/(2.0D0*dx%f(0))

       pi%s(0)  = (der_pi+der_psi)/2.0D0
       sigma%s(0) = (der_pi+der_psi)/2.0D0

       der_pi  = ( pi%f(grid_points-2)-4.0D0* pi%f(grid_points-1)+3.0D0 *pi%f(grid_points))/(2.0D0*dx%f(grid_points))
       der_psi = (sigma%f(grid_points-2)-4.0D0*sigma%f(grid_points-1)+3.0D0*sigma%f(grid_points))/(2.0D0*dx%f(grid_points))

       pi%s(grid_points) = -(der_pi-der_psi)/2.0D0
       sigma%s(grid_points) = (der_pi-der_psi)/2.0D0

    else

       print*,'Boundary condition not implemented!'
       stop
       
    endif

  end subroutine boundaries
