
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
  !   pi  - 3.1415......
  real(dp), parameter :: pi = 4.d0*datan(1.d0)

  ! floating point values
  !   time   - current time
  !   dt     - time step
  !   dt_max - maximal time step
  !   tend   - maximal time
  !   t_wo   - write out time
  !   cfl    - cfl number
  !   vel    - fuild velocity
  real(dp) :: time, dt, dt_max, tend, cfl 

end module globals




!******************************************************************************
! MODULE ADVECTION_HELPER: io, initialization & dt routines
!******************************************************************************
module convergence_helpers
  use globals
  use advection_io

  implicit none
  private

  ! public routines
  public :: read_param    ! read_param( filename{char} 
                          !   [, #headerlines{int}, skip_first_char{logical}] )      
  
  public :: write_conv    ! write_conv( folder{char}, initial_state{array(dp)}
                          !   [, suffix{char} -> default: dat] )

  public :: show_progress ! show_progres()
  
  public :: initialize()  ! initilaize (TODO: dependencies!)
  public :: calc_dt()     ! calc_dt(TODO: dependencies!)

  contains

    subroutine intialize_param(PHYmin,PHYmax,)
     
      implicit none
     
       
     
    end subroutine intialize_param
     
      
    subroutine initialize(state, state_ini, state_old, flux, dx)                
      ! Task: 
      !   allocates arrays, sets variables/parameters, set initial state
      ! Input:
      !   state     - state vector          
      !   state_ini - initial state vector    
      !   state_old - old state vector        
      !   flux      - flux vector             
      !   dx        - grid size             
      ! Outout:
      !   state     - state vector           -> fill   
      !   state_ini - initial state vector   -> alloc  
      !   state_old - old state vector       -> alloc  
      !   flux      - flux vector            -> alloc  
      !   dx        - grid size              -> calc
      !   time      - current time           -> set 0.0
      !   t_wo      - write out time         -> calc  
      !   PHYmax    - max index phys. domain -> calc   
      !   PHYmin    - min index phys. domain -> calc
      implicit none
  
      ! calculate/set values/parameters and allocate arrays
      call initialize_param(!TODO:stuff in here)
     
      ! set initial condition 
      




    end subroutine initialize



    subroutine calc_dt()
    
      implicit none
    
       
    
    end subroutine calc_dt
   
end module convergence_helpers


!******************************************************************************
! MODULE ADVECTION_FLUX: flux calculation -> rhs
!******************************************************************************
module advection_flux
  
  implicit none
  private

  ! parameters
  public :: flux  

  contains

    ! methods
    

end module advection_flux



!******************************************************************************
! MODULE RK2: runge kutta second order 
!******************************************************************************
module RK2
  
  implicit none
  save

  ! parameters
  public :: rk2

  contains

    ! methods
    <methods>

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
  !   state_old - old state (needed for RK2)
  !   state_ini - initial state
  !   flux      - flux vector
  float(dp), dimension(:), allocatable :: state, state_old, state_ini, flux

  ! floating point values
  !   dx   - grid size in x-direction
  !   t_wo - write out time
  float(dp) :: dx, t_wo
  
  
  call getarg(1, folder)

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


  do i = 1, 7

    ! improve resolution
    if (i > 2) then
      ncells = ncells**2
    end if

    write(*,*) '======================================================' 
    write(*,*) '    amount of cells:       ', ncells
    write(*,*) '======================================================'

    call info( '  initializing')
    call initialize(state,      &
                    state_ini,  &
                    state_old,  &
                    flux,       &
                    time,       &
                    t_wo,       &
                    PHYmax,     &
                    PHYmin,     &
                    dx)         
    call end_info()
      
    call info( '  copying initial state')
    state_ini = state                 
    call end_info() 


    ! integration
    
    do while (time <= tend)
      call rk2(state, alpha=0.5d0)
      call show_progress()
      
      if ( modulo(time,t_wo) == 0.0d0 ) then
        call write_conv(folder, state_ini)
      end if  
    end do


    write(*,*) ''
    write(*,*) '  Done!'
    write(*,*) ''
  end do


  write(*,*) '======================================================'
  write(*,*) color(' CONVERGENCE TEST SUCCESSFULLY DONE!', 'green')
  write(*,*) '======================================================'
    
end program conv_test
