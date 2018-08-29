!=================================================================================
! Simple Example: Strang Splitting of Adv-Diff Equation:  
!           phi_t + (u*phi - D*phi_x)_x = 0 (1D)
!           phi(0) = phi_0                  (IC)
!
! IC can be chosen via "inittype"	: 1-Gaussian, 2-Wave Packet
! Reconstruction Types "recontype"	: 1-Gudonov,  2-PLM:
!													|
!         											+--> Limiter: 1-MC,2-SuperBee
!=================================================================================

program OperatorSplitting_Simple

	implicit none
	
	!=============================================================================
	! Parameters
	!=============================================================================

	! TODO: Read Parameters from file
	! amount of cells and amount of ghost cells
	integer, parameter :: ncells 	= 64
	integer, parameter :: ngcells 	= 3

	! domain size
	double precision, parameter :: xmin = 0.d0
	double precision, parameter :: xmax = 1.d0

	! advection velocity
	double precision, parameter :: u = 1.d0

	! courant number cn = (u * dt / dx)
	double precision, parameter :: cn = 0.8d0

	! initial condition type
	integer, parameter :: inittype = 1
	
	! maximum simulation time
	double precision, parameter :: tmax = 5.d0

	
	double precision, dimension(2*ngcells + ncells) :: grid
	double precision, dimension(2*ngcells + ncells) :: phi, phi_left, phi_right, flux
	double precision, dimension(2*ngcells + ncells) :: phi_init

	double precision :: dx
	double precision :: t, dt 

	!=============================================================================
	! Program
	!=============================================================================

	! setup grid and initial conditions
	call setup_grid(ncells, ngcells, xmin, xmax, dx, grid)
	call init(ncells, ngcells, inittype, grid, phi)
	
	phi_init = phi

	t = 0.d0
	
	! write out initial state
	call write_out(args...)

	do while (t < tmax)
		
		! handle ghost cells
		call update_gcells(ncells, ngcells, phi)
		
		! calculate new time step
		call calc_dt             	
		if (t + dt > tmax) then  
			dt = tmax - t        
		end if

		! reconstruct phi
		call reconstruct(ncells, ngcells, dx, dt, u, recontype, limiter, phi, 		&
																phi_left, phi_right)

		! strang splitting: 
		! 	-> solve advection equation for dt/2 				--> phi* 
		!	-> solve diffusion equation for dt using phi*		--> phi**
		!	-> solve advection equation for dt/2 using phi**	--> phi(t+dt)
		call solve_advection()
		


		! update time
		t = t + dt
	end do 	
end program


!=================================================================================  
! Subroutines
!=================================================================================

!-------------------------+
! output: writes out data |
!-------------------------+
subroutine output(ncells, ngcells, inittype, recontype, limiter, u, t, grid, phi, folder)
	
	implicit none

	integer, intent(in) :: ncells, ngcells
	integer, intent(in) :: inittype, recontype, limiter
	
	double precision, intent(in) :: u, t
	
	double precision, dimension(2*ngcells+ncells), intent(in) :: grid, phi

	character(len=40), optional, intent(in) :: folder
	
	character(len=50) :: pre_string
	character(len=4)  :: time_string
	character(len=16) :: recon, init, res

	integer :: i, imin, imax

	imin = ngcells + 1
	imax = ngcells + ncells

	! generate recon, init strings
	if (recontype == 1) then
		recon = "godunov"
	else if (recontype == 2) then 
		if (limiter == 1) then
			recon = "plm+mc"
		else if (limiter == 2) then
			recon = "plm+SBee"
		end if
	end if

	if (inittype == 1) then
		init = "gaussian"
	else if (inittype == 2) then
		init = "packet"
	end if

	! open output file
	write(time_string, '(f4.2)') t
	write(res, '(i8)') ncells
	
	if (present(folder)) then
		pre_string="AdvDiff"//trim(folder)
	else
		pre_string="AdvDiff"
	end if

	open(unit=20, filename=pre_string//trim(recon)//"-"//trim(init)//"-ncells="//trim(adjustl(res))//"-t="//time_string, status="unknown")

	do i = imin, imax
		write(20, *) grid(i), phi(i)
	end do 

end subroutine output



!------------------------------+
! setup_grid: creates the grid |
!------------------------------+
subroutine setup_grid(ncells, ngcells, xmin, xmax, dx, grid)

	implicit none
		
	! amount of cells / ghost cells
	integer, intent(in) :: ncells, ngcells
		
	! maximum / minimum x value, x-step, grid array 
	double precision, intent(in) :: xmin, xmax
	double precision, intent(inout) :: dx 
	double precision, dimension(2*ngcells + ncells), intent(inout) :: grid

	! dummy index
	integer :: i 

	! create the grid
	dx = (xmax - xmin)/dble(ncells)

	do i = 1, 2*ngcells + ncells
		x(i) = (i - 0.5d0)*dx + xmin
	end do

	return
end subroutine setup_grid

!----------------------------------+
! init: set the initial conditions |
!----------------------------------+
subroutine init(ncells, ngcells, inittype, grid, phi)

	implicit none
		
	! number of cells / ghost cells
	integer, intent(in) :: ncells, ngcells

	! initialization type: 	1 -> Gaussian, 
	! 					   	2 -> Wave Paket
	integer, intent(in) :: inittype

	! grid, physical quantity phi
	double precision, dimension(2*ngcells+ncells), intent(in)		:: grid
	double precision, dimension(2*ngcells+ncells), intent(inout) 	:: phi

	! defining pi
	double precision, parameter :: pi = 4.d0*datan(1.d0)
		
	integer :: i, imin, imax
	
	imin = ngcells + 1
	imax = ngcells + ncells 		

	! loop over all cells and set initial conditions at the cell centre

	
	! TODO: Check +- 0.5d0 generic?
	if (inittype == 1) then
	! gaussian
		do i = imin, imax
			phi(i) = exp(-(x(i)- 0.5d0)**2/0.1d0**2)
		end do

	elseif(inittype == 2) then
	! wave paket
		do i = imin, imax
			phi(i) = sin(16.d0*pi*x(i))*exp(-36.d0*(x(i)-0.5d0)**2)
		end do
	end if 
				
	return
end subroutine init


!------------------------------------+
! update_gcells: updates ghost cells |
!------------------------------------+
subroutine update_gcells(ncells, ngcells, phi)

	implicit none

	integer, intent(in) :: ncells, ngcells	
	double precision, intent(inout) :: phi

	integer :: i, imax
	imax = ngcells + ncells

	! left boundary
	do i = 1, ngcells 
		phi(i) = phi(ncells+i)
	end do

	! right boundary
	do i = imax+1, 2*ngcells+ncells  
		phi(i) = phi(i-imax+ng)
	end do

	return
end subroutine update_gcells

!----------------------------------- +
! calc_dt: computes the new timestep |	
!------------------------------------+	
subroutine calc_dt(dx, u, cn, dt)

	implicit none

	double precision, intent(in) :: dx
	double precision, intent(in) :: u
	double precision, intent(in) :: cn
	double precision, intent(inout) :: dt
	
	! due to stability constrains
	dt = cn*dx/abs(u)

end subroutine calc_dt


!---------------------------------------------+
! reconstruct: computes the interfaces states |
!---------------------------------------------+
subroutine reconstruct(ncells, ngcells, dx, dt, u, recontype, limiter, &
		                                			& phi, phi_left, phi_right)

	implicit none

	integer, intent(in) :: ncells, ngcells
	integer, intent(in) :: recontype, limiter

	double precision, intent(in) :: dx, dt, u
	
	double precision, dimension(2*ngcells+ncells), intent(in) :: phi
	double precision, dimension(2*ngcells+ncells), intent(inout):: phi_left, phi_right

	double precision, dimension(2*ngcells+ncells) :: slope

	integer :: i, imin, imax

	imin = ngcells + 1
	imax = ngcells + ncells

	if (recontype == 1) then
		
		! Godunov's method (piecewise constant)
		
		! for each interface, we want to construct the left and the right 
		! states. Here, i referes to the left edge of zone i 

		! interfaces imin to imax+1 affect the data in cells [imin, imax]	
		do i = imin, imax+1
			! the left state on the current interface comes from cell i-1
			phi_left(i) = phi(i-1) 

			! the right state on the current interface comes from cell i 
			phi_right(i) = phi(i)
		end do
	
	elseif (recontype == 2) then

		! PLM with MC limiter
		if (limiter == 1) then

			! interface states are found by Taylor expansionin time (trough dt/2)
			! and space (dx/2 towards the interface).

			! for each interface left and right states have be reconstructed.
			! Here, interface i refers to the left edge of zone i.
		
			do i = imin-1, imax+1

				slope(i) = 	minmod( 											&
								minmod(	2.d0*( phi(i)   - phi(i-1) )/dx,   		&
										2.d0*( phi(i+1) - phi(i)   )/dx ),		&
								0.5d0*( phi(i+1) - phi(i-1) )/dx   				&
							) 															   
			end do
		
		! PLM with SuperBee limiter
		elseif (limiter == 2) then
			
			double precision :: slope1, slope2

			do i = imin-1, imax+1
				slope1 = minmod( (u(i+1) - u(i))/dx  , 2.d0*(u(i) - u(i-1))/dx )
				slope2 = minmod( 2.d0*(u(i+1) - u(i))/dx ,  (u(i) - u(i-1))/dx )

				slope(i) = maxmod(slope1, slope2)
			end do

		end if

		! interfaces imin to imax+1 affect the data in cells [imin, imax]
		do i = imin, imax+1
		
 			! the left state on the current interface comes from cell i-1
			phi_left(i)  = phi(i-1) + 0.5d0*dx*(1.d0 - u*(dt/dx))*slope(i-1)

			! the right state on the current interface comes from cell i 
			phi_right(i) = phi(i)   + 0.5d0*dx*(1.d0 + u*(dt/dx))*slope(i) 
		end do
 
	end if  

	return
end subroutine reconstruct


!---------------------------------------------------------------+
! solve_advection: calculates the fluxes for the advective part |
!---------------------------------------------------------------+
subroutine solve_advection(ncells, ngcells, u, dx, dt, phi, phi_left, phi_right, flux)

	implicit none

	integer, intent(in) :: ncells, ngcells
	
	double precision, dimension(2*ngcells+ncells), intent(in) :: phi_left, phi_right, flux
	double precision, dimension(2*ngcells+ncells), intent(inout) :: phi
	
	double precision, intent(in) :: u, dx, dt
	
	integer :: i, imin, imax

	imin = ngcells + 1
	imax = ngcells + ncells

	! loop over all interfaces and calculate the linear advective flux.   	
	! the advection  velocity tells us which direction is upwind.

	if (u <= 0.d0) then 
		do i = imin, imax+1
			flux(i) = u*phi_left(i)
		end do

	else
		do i = imin, imax+1
			flux(i) = u*phi_right(i)
		end do

	end if

	call update_advection(ncells, ngcells, dx, dt, phi, flux)
	return
end subroutine solve_advection



!-------------------------------------------------------------------+
! update: conservatively update the solution tho the new time level |
!-------------------------------------------------------------------+
subroutine update_advection(ncells, ngcells, dx, dt, phi, flux)

	implicit none
	
	integer, intent(in) :: ncells, ngcells
	
	double precision, intent(in) :: dx, dt
	
	double precision, dimension(2*ngcells+ncells), intent(in) :: flux
	double precision, dimension(2*ngcells+ncells), intent(inout) :: phi

	integer :: i, imin, imax


	imin = ngcells + 1
	imax = ngcells + ncells

	do i = imin, imax
		phi(i) = phi(i) + (dt/dx)*(f(i) - f(i+1))
	end do  

	return
end subroutine update_advection









!=================================================================================
! Functions
!=================================================================================

!---------------------------------------------------------------------------------
! Limiter Functions         
!---------------------------------------------------------------------------------

function minmod(a,b)
	
	implicit none

	double precision, intent(in)  :: a,b
	double precision, intent(out) :: minmod

	if (abs(a) < abs(b) .and. a*b > 0.d0) then
		minmod = a
	else if (abs(b) < abs(a) .and. a*b > 0.d0) then
		minmod = b
	else
		minmid = 0.d0
	endif

	return minmod

end function minmod


function maxmod(a,b)
 
    implicit none

    double precision, intent(in)  :: a,b
    double precision, intent(out) :: maxmod

    if (abs(a) > abs(b) .and. a*b > 0.d0) then
        maxmod = a
    else if (abs(b) > abs(a) .and. a*b > 0.d0) then
        maxmod = b
    else
        maxmid = 0.d0
    endif
 
    return maxmod
 
end function maxmod



