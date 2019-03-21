!******************************************************************************
! MODULE ADVECTION_FLUX: flux calculation -> rhs
!******************************************************************************
module advection_flux
  use globals
  use basic_helpers, only: boundaries
  use advection_reconstruction

  implicit none
  private

  interface fluxes_adv_upwind
	module procedure fluxes_adv_upwind_2D
	module procedure fluxes_adv_upwind_3D
  end interface

  ! public methods
  public :: fluxes_adv_upwind


  contains

    subroutine fluxes_adv_upwind_2D(state, fluxes, current_dim)
      ! Task:
      !   calculates the fluxes through interfaces using upwind scheme.
      ! Input:
      !   state       - actual state -> calculate interface values.
      !   current_dim - determines flux component.
      ! Output:
      !   fluxes - array containing computed fluxes.

      implicit none

	  integer , intent(in) :: current_dim
	  
      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells), intent(in)    :: state
      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells), intent(inout) :: fluxes

	  integer :: offsx, offsy
      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells,2) :: interfaces
	  real(dp) :: velocity
	  
	  ! set velocity as velocity in 'current_dim' direction.
	  velocity = vel(current_dim)
	  
	  if (velocity == 0.0_dp) then
		fluxes(:,:) = 0.0_dp
		return
	  end if  
	  
      ! reconstruct interface values
      call recon_lin_van_leer(state, interfaces, current_dim)
      !call godunov(state, interfaces) 

      ! set offsets
	  if (current_dim == 1) then
		offsx = -1
		offsy = 0
	  else if (current_dim == 2) then
		offsx = 0
		offsy = -1
	  end if
	  
      ! fluxes(i) correspond to the left interface of cell i [F(i-1/2)]
      if (velocity >= 0) then
        ! only right interface value matters
        fluxes(PHYmin:PHYmax, PHYmin:PHYmax) = velocity * &
			interfaces(PHYmin+offsx:PHYmax+offsx,PHYmin+offsy:PHYmax+offsy,2)
      else
        ! only left interface value matters
        fluxes(:,:) = velocity * interfaces(:,:,1)
      end if
      
	  call boundaries(fluxes)

    end subroutine fluxes_adv_upwind_2D



    subroutine fluxes_adv_upwind_3D(state, fluxes, current_dim)
      ! Task:
      !   calculates the fluxes through interfaces using upwind scheme.
      ! Input:
      !   state       - actual state -> calculate interface values.
      !   current_dim - determines flux component.
      ! Output:
      !   fluxes - array containing computed fluxes.

      implicit none

	  integer , intent(in) :: current_dim
	  
      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells,2*ngcells+ncells), intent(in) &
		:: state
      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells,2*ngcells+ncells), intent(inout)&
		:: fluxes

	  integer :: offsx, offsy, offsz
      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells,2*ngcells+ncells,2) :: interfaces
	  real(dp) :: velocity
	  
	  ! set velocity as velocity in 'current_dim' direction.
	  velocity = vel(current_dim)
	  
	  if (velocity == 0.0_dp) then
	    fluxes(:,:,:) = 0.0_dp
	    return
	  end if
	  
      ! reconstruct interface values
      call recon_lin_van_leer(state, interfaces, current_dim)
      !call godunov(state, interfaces) 

	  ! set offsets
	  if (current_dim == 1) then
		offsx = -1
		offsy = 0
		offsz = 0
	  else if (current_dim == 2) then
		offsx = 0
		offsy = -1
		offsz = 0
	  else if (current_dim == 3) then
		offsx = 0
		offsy = 0
		offsz = -1
	  end if
	  
      ! fluxes(i) correspond to the left interface of cell i [F(i-1/2)]
      if (velocity >= 0) then
        ! only right interface value matters
        fluxes(PHYmin:PHYmax,PHYmin:PHYmax, PHYmin:PHYmax) = velocity * &
			interfaces(PHYmin+offsx:PHYmax+offsx,PHYmin+offsy:PHYmax+offsy,PHYmin+offsz:PHYmax+offsz,2)
      else
        ! only left interface value matters
        fluxes(:,:,:) = velocity * interfaces(:,:,:,1)
      end if
	
	  call boundaries(fluxes)

    end subroutine fluxes_adv_upwind_3D

end module advection_flux
