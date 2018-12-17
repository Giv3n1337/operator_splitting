!******************************************************************************
! MODULE RECONSTRUCTION: reconstruct interface values
!******************************************************************************
module advection_reconstruction
  use globals
  use basic_helpers, only: boundaries

  implicit none
  save 

  private

  ! public methods
  public :: godunov
  public :: recon_lin_van_leer  ! paper ZIE04

  contains

    ! methods
    subroutine godunov(state, interfaces)
      ! Task:
      !   cell state is constant.
      ! Input:
      !   state      - actual state,
      !   interfaces - empty interfaces array. [dimension(:,2)]
      ! Output:
      !   interfaces - array containing left/right interface values.
      !                [interfaces(:,1) ->  left: state_{i-1/2}]
      !                [interfaces(:,2) -> right: state_{i+1/2}]

      implicit none

      real(dp), dimension(2*ngcells+ncells),   intent(in)    :: state
      real(dp), dimension(2*ngcells+ncells,2), intent(inout) :: interfaces
   
      ! dummy
      integer :: i

      do i = 1,2
        interfaces(:,i) = state(:)
      end do

      call boundaries(interfaces)
    end subroutine godunov
    

    subroutine recon_lin_van_leer(state, interfaces)
      ! Task:
      !   linear reconstruc left/right interfaces using van leer limiter.
      ! Input:
      !   state      - actual state,
      !   interfaces - empty interfaces array. [dimension(:,2)]
      !   ds_i       - delta_state_van_leer for cell i
      ! Output:
      !   interfaces - array containing reconstructed interfaces.
      !                [interfaces(:,1) ->  left: state_{i-1/2}]
      !                [interfaces(:,2) -> right: state_{i+1/2}]
  
      implicit none

      real(dp), dimension(2*ngcells+ncells),   intent(in)    :: state
      real(dp), dimension(2*ngcells+ncells,2), intent(inout) :: interfaces
    
      real(dp) :: ds

      ! dummy variable
      integer :: i

      ! fill interface values
      do i = PHYmin, PHYmax
          
          ! calculate ds
          ds = delta_state_van_leer( state(i-1:i+1) )    

          ! left/west interface
          interfaces(i,1) = state(i) - ds
          ! right/east interface
          interfaces(i,2) = state(i) + ds
      end do

      call boundaries(interfaces)
    end subroutine recon_lin_van_leer



    ! internal function calculating delta_state using van leer limiter
    function delta_state_van_leer(stencil) result(ds)

      implicit none

      real(dp), dimension(3), intent(in) :: stencil ! state(i-1:i+1)
      real(dp) :: ds

      ! van leer TVD limiter
      if ( abs(stencil(3) - stencil(1)) > 1.0d-16 ) then
        ds = max( (stencil(3)-stencil(2)) * (stencil(2)-stencil(1)), 0.0_dp )
        ds = ds / (stencil(3) - stencil(1))
      else 
        ds = 0.0_dp
      end if
    end function delta_state_van_leer

end module advection_reconstruction
