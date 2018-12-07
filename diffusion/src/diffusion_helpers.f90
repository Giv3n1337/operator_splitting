!******************************************************************************
! MODULE DIFFUSION_HELPERS: io, initialization & dt routines
!******************************************************************************
module diffusion_helpers
  use diffusion_io

  implicit none
  save

  private

  ! public routines
  public :: read_param      ! read_param( filename{char}
                            !   [, #headerlines{int}, skip_first_char{logical}] )

  public :: write_conv      ! write_conv( folder{char}, initial_state{array(dp)}
                            !   [, suffix{char} -> default: dat] )

  public :: show_progress   ! show_progress()
  public :: start_progress  ! start_progress()
  public :: end_progress    ! end_progress()
  public :: initialize      ! initilaize( state, state_ini, dt)


  contains

    subroutine initialize(state, state_ref, A_inv, dt)
      ! Task:
      !   allocates arrays, sets variables/parameters, set initial state
      ! Input:
      !   state     - state vector
      !   state_ref - reference state vector
      !   A_inv 	- inverse Matrix of A needed for crank nicolson scheme
      ! Outout:
      !   state     - state vector            -> fill
      !   state_ini - initial state vector    -> alloc
      !   dx        - grid size               -> calc
      !   dt        - time steps              -> set dt_max
      !   time      - current time            -> set 0.0
      !   t_wo      - write out time          -> calc
      !   PHYmax    - max index phys. domain  -> calc
      !   PHYmin    - min index phys. domain  -> calc
      !   A_inv     - inverse matrix of A depending on alpha
      
      implicit none

      ! floating point fields
      real(dp), dimension(:),   allocatable, intent(inout) :: state, state_ref
	  real(dp), dimension(:,:), allocatable, intent(inout) :: A_inv

      ! floating point values
      real(dp), intent(inout) :: dt

      ! dummy variable
      integer :: i

      ! calculate/set values/parameters and allocate arrays
      call initialize_param(state, state_ref, A_inv, dt)

      ! set initial condition
      do i = 1, ncells
        state(ngcells+i) = 1.0_dp / (sig*sqrt(2.0_dp*pi)) * &        ! Gaussian
          exp( -0.5_dp*( (i*dx - 0.5_dp*dx - 0.5_dp) / sig)**2 ) 

        ! if (i <= 0.25d0*ncells .or. i >= 0.75d0*ncells) then       ! Square
        !   state(ngcells+i) = 0.0d0
        ! else
        !   state(ngcells+i) = 2.0d0
        ! end if
      end do

      ! initialize progess bar
      call bar%initialize(                  &
        prefix_string='   progress ',       &
        bracket_left_string='|',            &
        bracket_left_color_fg='blue',       &
        bracket_right_string='| ',          &
        bracket_right_color_fg='blue',      &
        filled_char_string='+',             &
        filled_char_color_fg='yellow',      &
        empty_char_string=' ',              &
        add_progress_percent=.true.,        &
        progress_percent_color_fg='red',    &
        max_value=tend+0.0503_dp)          ! <-TODO:find better way

    end subroutine initialize


    ! internal methods

    subroutine initialize_param(state, state_ref, A_inv, dt)
      ! Task:
      !   allocates arrays, sets variables/parameters
      ! Input:
      !   state     - state vector
      !   state_ref - reference state vector
      !   A_inv 	- inverse matrix of A needed for crank nicolson
      !   dt        - time steps
      ! Outout:
      !   state     - state vector           -> alloc
      !   state_ref - reference state vector -> alloc (if not already)
      !   dx        - grid size              -> calc
      !   time      - current time           -> set 0.0
      !   t_wo      - write out time         -> calc
      !   dt        - time steps             -> set dt_max
      !   PHYmax    - max index phys. domain -> calc
      !   PHYmin    - min index phys. domain -> calc
      !   A_inv     - inverse matrix of A depending on alpha
      implicit none

      ! floating point fields
      real(dp), dimension(:),   allocatable, intent(inout) :: state, state_ref
      real(dp), dimension(:,:), allocatable, intent(inout) :: A_inv

      ! floating point values
      real(dp), intent(inout) :: dt

	  ! dummy
	  integer ::i
	  
	  ! info
	  integer :: info
	  
	  external DPOTR	! matrix inversion using cholesky decomposition
	  
      ! allocate arrays
      allocate(  state( ncells + 2*ngcells )  )
      allocate(  A_inv(ncells, ncells))
      
      if (.not. allocated(state_ref)) then
		allocate(  state_ref( ncells + 2*ngcells )  ) 
      end if

      ! calculate parameters
      dx     = 1.0d0 / real(ncells,dp)
      t_wo   =  tend !1.0d0 / vel
      PHYmin = ngcells + 1
      PHYmax = ngcells + ncells

      ! set time & step size
      time = 0.0_dp
      dt   = dt_max
      
      ! calculate alpha
      alpha = dt * diff_const / (2.0_dp*dx**2)
      
      ! construct A_inv
	  A_inv(:,:) = 0.0_dp
      do i = 1, ncells-1
		A_inv(i,i) = 1.0_dp + 2.0_dp * alpha
		A_inv(i,i+1) = -1.0_dp * alpha
		A_inv(i+1,i) = -1.0_dp * alpha
      end do
      
      A_inv(ncells,ncells) = 1.0_dp + 2.0_dp * alpha
	  A_inv(1, ncells)     = -1.0_dp * alpha
	  A_inv(ncells, 1)     = -1.0_dp * alpha
	  
	  call dpotri('U', ncells, A_inv, ncells, info)
 
    end subroutine initialize_param

end module diffusion_helpers
