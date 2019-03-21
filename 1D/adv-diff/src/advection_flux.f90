!******************************************************************************
! MODULE ADVECTION_FLUX: flux calculation -> rhs
!******************************************************************************
module advection_flux
  use globals
  use basic_helpers, only: boundaries
  use advection_reconstruction

  implicit none
  private

  ! public methods
  public :: fluxes_adv_upwind


  contains

    subroutine fluxes_adv_upwind_2D(state, fluxes, current_dim)
      ! Task:
      !   calculates the fluxes through interfaces using upwind scheme.
      ! Input:
      !   state - actual state -> calculate interface values,
      !   current_dim - current dimension.
      ! Output:
      !   fluxes - array containing computed fluxes.

      implicit none

      real(dp), dimension(:,:), intent(in) :: state
      
      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells), &
        intent(inout) :: fluxes

      real(dp), dimension(2*ngcells+ncells,2) :: interfaces

      ! reconstruct interface values
      call recon_lin_van_leer(state, interfaces)
      !call godunov(state, interfaces) 

      ! fluxes(i) correspond to the left interface of cell i [F(i-1/2)]
      if (vel >= 0) then
        ! only right interface value matters
        fluxes(PHYmin:PHYmax) = vel * interfaces(PHYmin-1:PHYmax-1,2)
        call boundaries(fluxes)
      else
        ! only left interface value matters
        fluxes(:) = vel * interfaces(:,1)
      end if

    end subroutine fluxes_adv_upwind


end module advection_flux
