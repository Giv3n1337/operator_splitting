!******************************************************************************
! MODULE RECONSTRUCTION: reconstruct interface values
!******************************************************************************
module advection_reconstruction
  use globals
  use basic_helpers, only: boundaries

  implicit none
  save 

  private

  interface recon_lin_van_leer
    module procedure recon_lin_van_leer_2D
    module procedure recon_lin_van_leer_3D
  end interface

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
    
    
 
    subroutine recon_lin_van_leer_2D(state, interfaces, current_dim)
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

	  integer, intent(in) :: current_dim

      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells),   intent(in)    &
		:: state
      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells,2), intent(inout) &
		:: interfaces
		
		
	  integer :: offsx, offsy
      real(dp) :: ds

      ! dummy variable
      integer :: i, j

      ! set offsets
      if (current_dim == 1) then
		offsx = 1
		offsy = 0
	  else if (current_dim == 2) then
		offsx = 0
		offsy = 1
	  end if

      ! fill interface values
      do j = PHYmin, PHYmax
        do i = PHYmin, PHYmax  
          ! calculate ds
          ds = delta_state_van_leer( state(i-offsx:i+offsx,j-offsy:j+offsy) )    

          ! left/west interface
          interfaces(i,j,1) = state(i,j) - ds
          ! right/east interface
          interfaces(i,j,2) = state(i,j) + ds
        end do
      end do

      call boundaries(interfaces)
      
    end subroutine recon_lin_van_leer_2D



    subroutine recon_lin_van_leer_3D(state, interfaces, current_dim)
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

	  integer, intent(in) :: current_dim

      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells,2*ngcells+ncells)  ,&
		intent(in) :: state
      real(dp), dimension(2*ngcells+ncells,2*ngcells+ncells,2*ngcells+ncells,2),&
        intent(inout) :: interfaces
		
	  integer :: offsx, offsy, offsz
      real(dp) :: ds

      ! dummy variable
      integer :: i, j, k

      ! set offsets
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

      ! fill interface values
      do k = PHYmin, PHYmax
        do j = PHYmin, PHYmax
          do i = PHYmin, PHYmax
            
			! calculate ds
			ds = delta_state_van_leer( &
			  state(i-offsx:i+offsx, j-offsy:j+offsy, k-offsz:k+offsz) )    

			! left/west interface
			interfaces(i,j,k,1) = state(i,j,k) - ds
			! right/east interface
			interfaces(i,j,k,2) = state(i,j,k) + ds
		  end do
        end do
      end do

      call boundaries(interfaces)
    end subroutine recon_lin_van_leer_3D
    

    ! internal function calculating delta_state using van leer limiter
    function delta_state_van_leer(stencil) result(ds)

      implicit none

      real(dp), dimension(3), intent(in) :: stencil 
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
