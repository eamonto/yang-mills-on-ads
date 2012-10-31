
! ===========================================================================
! main.f90
! ===========================================================================
! Principal program, this program solves the Yang-Mills equation in 1+1 dimensions.

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


  program main

    use param
    use mylibrary

    implicit none

    integer :: l,k

    real(double) :: time = initial_time

    type(dynamical_func) omega,pi,sigma
    type(extra_func)     x,dx   
    type(scalar_func)    energy

    !Allocation of function that are evolve
    call allocate_dyn(omega,Nx)  !NO USER
    call allocate_dyn(pi   ,Nx)  !NO USER
    call allocate_dyn(sigma,Nx)  !NO USER

    !Allocation of function that are not evolve
    call allocate_extra(x ,Nx)  !NO USER
    call allocate_extra(dx,Nx)  !NO USER

    !Initialization of the functions, including the grid
    call initial_data(omega,pi,sigma,x,dx,Nx)  !USER

    !Compute the energy
    call compute_observables(energy,omega,pi,sigma,x,dx,Nx) !USER

    !Created the output file for the dynamical function in  
    !the "output_dir" with name "output_file" and id "file_number"
    call create_output_dyn(omega,output_dir,"omega.x",200)     !NO USER
    call create_output_dyn(pi   ,output_dir,"pi.x"   ,201)     !NO USER
    call create_output_dyn(sigma,output_dir,"sigma.x",202)     !NO USER

    !Created the output file for the extra function in  
    !the "output_dir" with name "output_file" and id "file_number"
    call create_output_extra(dx,output_dir,"dx.x",203)     !NO USER

    !Created the output file for the scalar function in  
    !the "output_dir" with name "output_file" and id "file_number"
    call create_output_scalar(energy,output_dir,"energy.t",204)  !NO USER

    !Print the output of two vector functions on the grid
    call output_obs_obs(x%f,omega%f,omega%name,omega%id,Nx) !NO USER
    call output_obs_obs(x%f,pi%f   ,pi%name   ,pi%id   ,Nx) !NO USER
    call output_obs_obs(x%f,sigma%f,sigma%name,sigma%id,Nx) !NO USER

    call output_obs_obs(x%f,dx%f,dx%name,dx%id,Nx) !NO USER

    call output_obs(time,energy) !NO USER
 
    print *,''
    print *,'    Time     '
    print *,'-------------'
    write(*,"(F12.5)") time

    do l=1,Ntime

       !Store the initial values before the integration
       call store_levels_rk4(omega)  !NO USER
       call store_levels_rk4(pi)     !NO USER
       call store_levels_rk4(sigma)  !NO USER

       do k=1,4

          !Derivatives of the functions to be evolve
          call sources(omega,pi,sigma,x,dx,Nx)               !USER

          !Boundaries
          call boundaries(omega,pi,sigma,x,dx,Nx,boundary)   !USER

          !Implemetation of the Runge-Kutta 4 method 
          call evolution_rk4(k,omega,dt)  !NO USER
          call evolution_rk4(k,pi   ,dt)  !NO USER
          call evolution_rk4(k,sigma,dt)  !NO USER

       enddo

       time = time + dt

       if(mod(l,every_0D).eq.0) then
          write(*,"(F12.5)") time

          !Compute the energy
          call compute_observables(energy,omega,pi,sigma,x,dx,Nx) !USER

          call output_obs(time,energy) !NO USER
       endif

       !Print the output of two scalar functions on the grid
       if(mod(l,every_1D).eq.0) then
          call output_obs_obs(x%f,omega%f,omega%name,omega%id,Nx)  !NO USER
          call output_obs_obs(x%f,pi%f   ,pi%name   ,pi%id   ,Nx)  !NO USER
          call output_obs_obs(x%f,sigma%f,sigma%name,sigma%id,Nx)  !NO USER
       endif

    enddo

    !Deallocation of function that are evolve
    call deallocate_dyn(omega)  !NO USER
    call deallocate_dyn(pi)     !NO USER
    call deallocate_dyn(sigma)  !NO USER

    !Deallocation of function that are not evolve
    call deallocate_extra(x)   !NO USER
    call deallocate_extra(dx)  !NO USER

    print *,''
    print *,'   Finish'
    print *,''

  end program main
