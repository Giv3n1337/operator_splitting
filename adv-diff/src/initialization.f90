!******************************************************************************
! MODULE DIFFUSION_HELPERS: io, initialization & dt routines
!******************************************************************************
module initialization
  use convergence_test_io

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

    subroutine initialize(state, state_ref, dt)
      ! Task:
      !   allocates arrays, sets variables/parameters, set initial state
      ! Input:
      !   state     - state vector
      !   state_ref - reference state vector.
      ! Output:
      !   state     - state vector            -> fill
      !   state_ini - initial state vector    -> alloc
      !   dx        - grid size               -> calc
      !   dt        - time steps              -> set dt_max
      !   time      - current time            -> set 0.0
      !   t_wo      - write out time          -> calc
      !   PHYmax    - max index phys. domain  -> calc
      !   PHYmin    - min index phys. domain  -> calc 
      !   A_inv     - inverse coeff matrix A  -> alloc   |_ needed for 
      !   B         - coeff matrix B          -> alloc   |  crank-nicolson
      implicit none

      ! floating point fields
      real(dp), dimension(:),   allocatable, intent(inout) :: state, state_ref

      ! floating point values
      real(dp), intent(inout) :: dt

      ! dummy variable
      integer :: i

      ! calculate/set values/parameters and allocate arrays
      call initialize_param(state, state_ref, dt)

      ! set initial condition and reference solution
      do i = 1, ncells
        !state(ngcells+i) = 1.0_dp / (sig*sqrt(2.0_dp*pi)) * &      ! Gaussian
        !  exp( -0.5_dp*( (i*dx - 0.5_dp*dx - 0.5_dp) / sig)**2 ) 

        ! if (i <= 0.25d0*ncells .or. i >= 0.75d0*ncells) then      ! Square
        !   state(ngcells+i) = 0.0d0
        ! else
        !   state(ngcells+i) = 2.0d0
        ! end if
        
        state(ngcells+i) = sin(2.0_dp*pi*(i-0.5_dp)*dx) 		      	! Sin
        
        
        state_ref(ngcells+i) = sin(2.0_dp*pi*(i-0.5_dp)*dx) *    &  ! ref
			    exp(-1.0_dp * diff_const* (2.0_dp * pi)**2 * tend)
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

    subroutine initialize_param(state, state_ref, dt)
      ! Task:
      !   allocates arrays, sets variables/parameters
      ! Input:
      !   state     - state vector,
      !   state_ref - reference state vector,
      !   dt        - time steps.
      ! Outout:
      !   state     - state vector           -> alloc
      !   state_ref - reference state vector -> alloc (if not already)
      !   dx        - grid size              -> calc
      !   time      - current time           -> set 0.0
      !   t_wo      - write out time         -> calc
      !   dt        - time steps             -> set dt_max
      !   PHYmax    - max index phys. domain -> calc
      !   PHYmin    - min index phys. domain -> calc
      !   A_inv     - inverse coeff matrix A -> alloc   |_ needed for 
      !   B         - coeff matrix B         -> alloc   |  crank-nicolson
      implicit none

      ! floating point fields
      real(dp), dimension(:),   allocatable, intent(inout) :: state, state_ref
      
      ! floating point values
      real(dp), intent(inout) :: dt

	    ! dummy
	    integer ::i
	  
      ! allocate arrays
      allocate(      state(ncells+2*ngcells )  )
      allocate(  state_ref(ncells+2*ngcells )  ) 
      allocate(      A_inv(ncells, ncells)  )
      allocate(          B(ncells, ncells)  )
      
      ! calculate parameters
      dx     = 1.0d0 / real(ncells,dp)
      t_wo   =  tend !1.0d0 / vel
      PHYmin = ngcells + 1
      PHYmax = ngcells + ncells

      ! set time & step size
      time = 0.0_dp
      dt   = dt_max
      
    end subroutine initialize_param

end module initialization
