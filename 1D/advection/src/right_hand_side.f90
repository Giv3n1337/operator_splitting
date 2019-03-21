!******************************************************************************
! MODULE RHS: generates rhs using the calculated fluxes
!******************************************************************************
module right_hand_side
  use globals
  use flux

  implicit none
  private

  ! public methods
  public :: rhs

  contains

    subroutine rhs(state, res)
      ! Task:
      !   calculates the right hand side d/dt(state).
      ! Input:
      !   state - actuall state vector,
      !   res   - resulting rhs vector (empty).
      ! Output:
      !   res   - calculated rhs vector.

      implicit none

      real(dp), dimension(ncells+2*ngcells), intent(in)    :: state
      real(dp), dimension(ncells),           intent(inout) :: res

      real(dp), dimension(ncells+2*ngcells) :: fluxes

      ! compute fluxes through interfaces
      call fluxes_adv_upwind(state, fluxes)

      ! calculating rhs
      res(:) = (fluxes(PHYmin:PHYmax) - fluxes(PHYmin+1:PHYmax+1)) / dx
    end subroutine rhs

end module right_hand_side
