program main
  use advection

  implicit none

  character(len=:), allocatable :: folder
 
  call getarg(1, folder)

  ! set up simulation
  write(*,*) 'Setting up simulation:'
  write(*,*) '======================================================'
  write(*,*) ''
  call info("  reading in parameters from 'parameters.txt'")
  call read_param('parameters.txt', skip_first_char==.true.)
  call end_info()
  write(*,*) '    amount of cells:  ', ncells
  write(*,*) '    maximal time:     ', time_end
  write(*,*) '    maximal timestep: ', dt_max 
  write(*,*) '    cfl number:       ', cfl

  call info('  initializing')
  call initialize(time, state, dx, flux, prop_speed, state_tmp, state_old)
  call end_info()

  call info('  write out initial state')
  call write_state(folder)
  call end_info()
  write(*,*) ''

  write(*,*) 'Starting integration:'
  write(*,*) '======================================================' 
  do while (time <= time_end)    
    call rk2(state, alpha=0.5d0)
    call show_progress()
    call write_state(folder)  
  end do

  write(*,*) ''
  write(*,*) color('SIMULATION SUCCESSFULLY DONE!', 'green')
  write(*,*) ''
end program main
