!******************************************************************************
! MODULE ADVECTION_RHS: generates rhs using the calculated fluxes
!******************************************************************************
module advection_rhs
  use globals
  use advection_flux

  implicit none
  private

  interface rhs
    module procedure rhs_2D
    module procedure rhs_3D
  end interface

  ! public methods
  public :: rhs

  contains

    subroutine rhs_2D(s, ds_dt)
      ! Task:
      !   calculates the right hand side d/dt(state).
      ! Input:
      !   s     - actuall state vector,
      !   ds_dt - resulting rhs vector (empty).
      ! Output:
      !   ds_dt - calculated rhs.

      implicit none

      real(dp), dimension(:,:), intent(in)    :: s
      real(dp), dimension(:,:), intent(inout) :: ds_dt

      real(dp), dimension(ncells+2*ngcells,ncells+2*ngcells) :: fluxes
  
      integer :: current_dim, offsx, offsy

      ! dummy
      integer :: i,j

      ds_dt(:,:) = 0.0_dp

      do current_dim = 1,2 
        ! compute fluxes through interfaces
        call fluxes_adv_upwind(s, fluxes, current_dim)

        if (current_dim == 1) then
          offsx = 1 
          offsy = 0
        else if (current_dim == 2) then
          offsx = 0
          offsy = 1
        end if

        ! calculating rhs
        do j = PHYmin, PHYmax
          do i = PHYmin, PHYmax
            ds_dt(i,j) = ds_dt(i,j) + & 
              (fluxes(i,j) - fluxes(i+offsx,j+offsy)) / dx
          end do
        end do
        
      end do
    end subroutine rhs_2D

    
    subroutine rhs_3D(s, ds_dt)
      ! Task:
      !   calculates the right hand side d/dt(state).
      ! Input:
      !   s     - actuall state vector,
      !   ds_dt - resulting rhs vector (empty).
      ! Output:
      !   ds_dt - calculated rhs.

      implicit none

      real(dp), dimension(:,:,:), intent(in)    :: s
      real(dp), dimension(:,:,:), intent(inout) :: ds_dt

      real(dp), dimension(ncells+2*ngcells,ncells+2*ngcells,ncells+2*ngcells) &
		:: fluxes
  
      integer :: current_dim, offsx, offsy, offsz

      ! dummy
      integer :: i,j,k

      ds_dt(:,:,:) = 0.0_dp

      do current_dim = 1,3
        ! compute fluxes through interfaces
        call fluxes_adv_upwind(s, fluxes, current_dim)

        if (current_dim == 1) then
          offsx = 1 
          offsy = 0
          offsz = 0
        else if (current_dim == 2) then
          offsx = 0
          offsy = 1
          offsz = 0
        else if (current_dim == 3) then
          offsx = 0
          offsy = 0
          offsz = 1
        end if

        ! calculating rhs
        do k = PHYmin, PHYmax
          do j = PHYmin, PHYmax
            do i = PHYmin, PHYmax
              ds_dt(i,j,k) = ds_dt(i,j,k) + & 
                (fluxes(i,j,k) - fluxes(i+offsx,j+offsy,k+offsz)) / dx
            end do
          end do
        end do

      end do
    end subroutine rhs_3D

end module advection_rhs
