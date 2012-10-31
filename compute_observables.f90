
! ===========================================================================
! compute_observables.f90
! ===========================================================================
! Compute observables (function define by user, e.g. Energy )

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

! Up to date: 30 Oct 2012					


  subroutine compute_observables(energy,omega,pi,sigma,x,dx,grid_points)

    use param
    use mylibrary

    implicit none

    type(scalar_func) energy
    type(dynamical_func) omega,pi,sigma
    type(extra_func)     x,dx   

    integer :: grid_points

    integer :: i

    dx%f = dx_aux

    energy%f = 0.0D0

    do i=1,grid_points
       energy%f = energy%f + dx%f(i)*( pi%f(i)*pi%f(i) + sigma%f(i)*sigma%f(i) &
                     + (omega%f(i)*omega%f(i)-1.0D0)**2/(2.0D0*sin(x%f(i))**2) )
    enddo

  end subroutine compute_observables
