
! ===========================================================================
! param.f90
! ===========================================================================
! Global parameters for the physical system.

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


  module param

    use mylibrary

    implicit none

    character(100), parameter :: output_dir = "test"  !output directory

    integer, parameter :: boundary = 1                 !1=zero 2=eigenvalue

    integer, parameter :: Nx = 1000                    !Number of points in the grid

    real(double), parameter :: courant = 0.5D0         !Courant factor
    real(double), parameter :: dx_aux = 3.14159265359D0/(2.0D0*Nx) !dx homogeneous
    real(double), parameter :: dt = courant*dx_aux     !time step

    real(double), parameter :: initial_time =0.0D0     !initial time
    real(double), parameter :: total_time =1.5D0       !final time
    integer     , parameter :: Ntime = int(total_time/dt)  !number of time steps

    integer     , parameter :: every_0D = 1            !time output
    integer     , parameter :: every_1D = int(Ntime/10)!spatial output

  end module param
