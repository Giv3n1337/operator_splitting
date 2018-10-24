
!******************************************************************************
! MODULE GLOBALS: global parameters and variables
!******************************************************************************
module globals 

  implicit none
  save  

  ! intege parameters
  !   dp - double precision
  integer, parameter :: dp = selected_real_kind(15,307)

  ! integer values
  !   ncells  - number of cells
  !   ngcells - number of ghost cells
  !   PHYmin  - first index of physical domain
  !   PHYmax  - last index of physical domain
  integer :: ncells, ngcells, PHYmin, PHYmax

  ! floating point parameters
  !   pi - 3.1415......
  real(dp), parameter :: pi = 4.d0*datan(1.d0)

  ! floating point values
  !   time   - current time
  !   dt_max - maximal time step
  !   tend   - maximal time
  !   t_wo   - write out time
  !   cfl    - cfl number
  !   vel    - fuild velocity
  !   dx     - cell size
  real(dp) :: time, dt_max, tend, cfl, vel, dx, t_wo

end module globals




!******************************************************************************
! MODULE ADVECTION_HELPER: io, initialization & dt routines
!******************************************************************************
module convergence_helpers
  use advection_io

  implicit none
  private

  ! public routines
  public :: read_param    ! read_param( filename{char} 
                          !   [, #headerlines{int}, skip_first_char{logical}] )      
  
  public :: write_conv    ! write_conv( folder{char}, initial_state{array(dp)}
                          !   [, suffix{char} -> default: dat] )

  public :: show_progress ! show_progres()


  public :: initialize    ! initilaiz( state, state_ini, dt)


  contains 
      
    subroutine initialize(state, state_ini, dt)                
      ! Task: 
      !   allocates arrays, sets variables/parameters, set initial state
      ! Input:
      !   state     - state vector          
      !   state_ini - initial state vector    
      ! Outout:
      !   state     - state vector           -> fill   
      !   state_ini - initial state vector   -> alloc  
      !   dx        - grid size              -> calc
      !   dt        - time steps             -> set dt_max
      !   time      - current time           -> set 0.0
      !   t_wo      - write out time         -> calc  
      !   PHYmax    - max index phys. domain -> calc   
      !   PHYmin    - min index phys. domain -> calc
      
      implicit none
  
      ! floating point fields
      real(dp), dimension(:), allocatable, intent(inout) :: state, state_ini
      
      ! floating point values
      real(dp), intent(inout) :: dx, t_wo, dt

      ! calculate/set values/parameters and allocate arrays
      call _initialize_param(state, state_ini, t_wo, dt, dx)
     
      ! set initial condition 
      do i = PHYmin, PHYmax
        state(i) = 1.d0 / (0.25d0*sqrt(2.0d0*pi)) * &  ! Gaussian
          exp( -0.5d0*( (i - 0.5d0) / (0.25d0) )**2

        ! if (i <= 0.25d0 .or. i >= 0.75d0) then       ! Square
        !   state(i) = 0.0_dp
        ! else
        !   state(i) = 2.0d0
        ! end if
      end do

      ! initialize progess bar
      call bar%initialize(
        prefix_string='  progress ', &
        bracket_left_string='[',     &
        bracket_right_string=']',    &
        filled_char_string='=',      &
        empty_char_string=' ',       &
        add_progress_procent=.true.)
      
    end subroutine initialize


    ! internal methods

    subroutine _initialize_param(state, state_ini, dt)
      ! Task: 
      !   allocates arrays, sets variables/parameters
      ! Input:
      !   state     - state vector          
      !   state_ini - initial state vector                         
      !   dt        - time steps      
      ! Outout:
      !   state     - state vector           -> alloc  
      !   state_ini - initial state vector   -> alloc  
      !   dx        - grid size              -> calc
      !   time      - current time           -> set 0.0
      !   t_wo      - write out time         -> calc  
      !   dt        - time steps             -> set dt_max
      !   PHYmax    - max index phys. domain -> calc   
      !   PHYmin    - min index phys. domain -> calc
      
      implicit none
     
      ! floating point fields
      real(dp), dimension(:), allocatable, intent(inout) :: state, state_ini
      
      ! floating point values
      real(dp), intent(inout) :: dt

      ! allocate arrays
      allocate(                           &
            state( ncells + 2*ngcells ),  &
        state_ini( ncells + 2*ngcells )   &
      )   

      ! calculate parameters
      dx     = 1.0d0 / real(ncells,dp)
      t_wo   = 1.0d0 / vel
      PHYmin = ngcells + 1
      PHYmax = ngcells + ncells

      ! set time & step size
      time = 0.0_dp
      dt   = dt_max

    end subroutine initialize_param

end module convergence_helpers




!******************************************************************************
! MODULE BASIC_HELPERS: contains basic methods andinterfaces
!******************************************************************************
module basic_helpers
  
  implicit none
  private

  ! pulbic methods
  !public :: calc_dt      <- TODO: not needed atm.
  public :: boundaries    ! boundaries(array)

  interface boundaries
    module procedure boundaries_1D
    module procedure boundaries_nD
  end interface

  contains 

    subroutine boundaries_1D(state)
      ! Task:
      !   updates periodic boundary conditions of a 1D array.
      ! Input:
      !   state - state vector.
      ! Output:
      !   state - updated state vector.

      implicit none
   
      real(dp), dimension(:), intent(inout) :: state
        
      state(1:ngcells) = state(PHYmax+1-ngcells:PHYmax)
      state(PHYmax+1:) = state(PHYmin:PHYmin-1+ngcells)
   
    end subroutine boundaries_1D
   

    subroutine boundaries_nD(state)
      ! Task:
      !   updates periodic boudary conditions of a n-dim array.
      ! Input:
      !   state - state vector
      ! Output:
      !   state - updated state vector.

      implicit none

      real(dp), dimension(:,:), intent(inout) :: state
      integer :: i
      
      do i=1, size(state,2)
        state(1:ngcells,i) = state(PHYmax+1-ngcells:PHYmax,i)
        state(PHYmax+1:,i) = state(PHYmin:PHYmin-1+ngcells,i)
      end do

    end subroutine boundaries_nD


end module basic_helpers




!******************************************************************************
! MODULE RECONSTRUCTION: reconstruct interface values
!******************************************************************************
module reconstruction
  use globals
  use basic_helpers, only: boundaries

  implicit none
  private

  ! public methods
  public :: recon_lin_van_leer(state, interfaces)  ! paper ZIE04

  contains
    
    ! methods
    subroutine recon_lin_van_leer(state, interfaces)
      ! Task:
      !   linear reconstruc left/right interfaces using van leer limiter.
      ! Input:
      !   state      - actual state,
      !   interfaces - empty interfaces array. [dimension(:,2)]
      ! Output:
      !   interfaces - array containing reconstructed interfaces.
      !                [interfaces(:,1) ->  left: state_{i-1/2}]
      !                [interfaces(:,2) -> right: state_{i+1/2}]
      implicit none

      real(dp), dimension(:),   intent(in)    = state
      real(dp), dimension(:,2), intent(inout) = interfaces = 0.0_dp


      real(dp) = delta_state

      ! fill interface values
      do i = PHYmin, PHYmax 
          delta_state = _delta_state_van_leer(state, i)

          ! left/west interface
          interfaces(i,1) = state(i) - delta_state_van_leer(state, i)
          ! right/east interface  
          interfaces(i,2) = state(i) + delta_state_van_leer(state, i)
      end do

      call boundaries(interfaces)  
    end subroutine recon_lin_van_leer



    ! internal function calculating delta_state using van leer limiter
    function _delta_state_van_leer(state, i) result(ds)

      implicit none 

      real(dp), dimension(:), intent(in) :: state ! state vector
      integer, intent(in) :: i                    ! actual index

      real(dp), intent(out) ds

      ! van leer TVD limiter
      ds = max( (state(i+1)-state(i)) * (state(i)-state(i-1) ), 0.0_dp )
      ds = ds / (state(i+1) - state(i-1))
    end function _delta_state_van_leer

end module reconstruction




!******************************************************************************
! MODULE ADVECTION_FLUX: flux calculation -> rhs
!******************************************************************************
module flux
  use reconstruction

  implicit none
  private

  ! public methods
  public :: fluxes_adv_upwind  


  contains

    subroutine fluxes_adv_upwind(state, fluxes)
      ! Task:
      !   calculates the fluxes through interfaces using upwind scheme.
      ! Input:
      !   state  - actual state -> calculate interface values.
      ! Output:
      !   fluxes - array containing computed fluxes.

      implicit none

      real(dp), dimension(:), intent(in)    :: state
      real(dp), dimension(:), intent(inout) :: fluxes

      real(dp), dimension(:,2) :: interfaces

      ! reconstruct interface values
      call recon_lin_van_leer(state, interfaces)

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
    

end module flux



!******************************************************************************
! MODULE RHS: generates rhs using the calculated fluxes
!******************************************************************************
module rhs
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

      real(dp), dimension(ncells+2*ngcell) :: fluxes

      ! compute fluxes through interfaces
      call fluxes_adv_upwind(state, fluxes)
      
      ! calculating rhs
      res(:) = (fluxes(PHYmin:PHYmax) - fluxes(PHYmin+1:PHYmax+1)) / dx

    end subroutine rhs 

end module rhs




!******************************************************************************
! MODULE RK2: runge kutta second order 
!******************************************************************************
module RungeKutta
  use rhs

  implicit none
  private


  public :: rk2
  
  contains

    subroutine rk2(state, dt, const_dt)
      ! Task:
      !   second order runge kutta method.
      ! Input:
      !   state    - actual state vector,
      !   dt       - last time step,
      !   const_dt - indicator, whether dt = const. [0: false, 1: true]
      ! Output:
      !   state - updated state vector,
      !   dt    - actual time step.    
     
      implicit none
    
      integer, intent(in) :: const_dt 

      real(dp), intent(inout) :: dt
      real(dp), dimension(ncells+2*ngcells), intent(inout) :: state

      ! runge kutta parameters
      real(dp), dimension(ncells+2*ngcells) :: k1
      
      ! rhs of d/dt(state) [first RK-step], rhs of d/dt(k1) [second RK-step]
      real(dp), dimension(ncells) :: k2

      ! update boundary conditions
      call boundaries(state)

      ! runge-kutta steps
      if (const_dt == 0) then
        call update_dt
        call rhs(state, k2, dt)
      else
        call rhs(state, k2, dt)
      end if  

      k1(PHYmin:PHYmax) = state(PHYmin:PHYmax) + dt*k2

      call boundaries(k1)
      call rhs(k1, k2, dt)

      state(PHYmin:PHYmax) = &
        0.5_dp * ( state(PHYmin:PHYmax) + k_1(PHYmin:PHYmax) + dt*k2 )

    end subroutine rk2
  

end module RK2





!******************************************************************************
! MAIN PROGRAM
!******************************************************************************
program conv_test
  use globals
  use convergence_helpers
  use RK2

  implicit none

  ! characters
  !   folder - where to store output files
  character(len=:), allocatable = folder 

  ! floating point arrays
  !   state     - state vector
  !   state_ini - initial state 
  float(dp), dimension(:), allocatable :: state, state_ini

  ! floating point values
  !   dx   - grid size in x-direction
  !   t_wo - write out time
  !   dt   - current time step
  float(dp) :: dx, t_wo, dt

  ! set folder 
  if iargc() = 1
    call getarg(1, folder)
  else 
    folder = fullname('.')
  end if

  !*******************
  ! CONVERGENCE TEST 
  !*******************

  ! initialization                                            
  write(*,*) '======================================================' 
  write(*,*) ' Setting up Convergence Test:'
  write(*,*) '======================================================'
  !write(*,*) ''
  write(*,*) '  output folder: ', folder

  call info( "  reading in parameters from 'parameters.txt'")  
  call read_param('parameters.txt', skip_first_char=.true.)
  call end_info()
  write(*,*) '    amount of cells:       ', ncells
  write(*,*) '    amount of ghost cells: ', ngcells
  write(*,*) '    maximal time:          ', tend
  write(*,*) '    maximal timestep:      ', dt_max
  write(*,*) '    cfl number:            ', cfl
  write(*,*) '    velocity:              ', vel


  ! loop over resplutions
  do i = 1, 7

    ! improve resolution
    if (i > 2) then
      ncells = ncells**2
    end if

    write(*,*) '======================================================' 
    write(*,*) '    amount of cells:       ', ncells
    write(*,*) '======================================================'

    call info( '  initializing')
    call initialize(state, state_ini, dt)         
    call end_info()
      
    call info( '  copying initial state')
    state_ini(:) = state(:)                 
    call end_info() 


    ! integration
    call start_progress()
    
    do while (time <= tend)
      call rk2(state, dt, 1)
      call show_progress()
      
      if ( modulo(time,t_wo) == 0.0d0 ) then
        call write_conv(folder, state_ini)
      end if
      
      time = time + dt  
    end do


    write(*,*) ''
    write(*,*) '  Done!'
    write(*,*) ''
  end do


  write(*,*) '======================================================'
  write(*,*) color(' CONVERGENCE TEST SUCCESSFULLY DONE!', 'green')
  write(*,*) '======================================================'
    
end program conv_test
