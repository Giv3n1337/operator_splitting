!******************************************************************************
! MAIN PROGRAM: convergence test of diffusion equation
!******************************************************************************
program conv_test
  use diffusion_io
  use diffusion_helpers
  use crank_nicolson
  
  implicit none

  ! characters
  !   folder - where to store output files,
  !   arg    - arguments passed in command line when executeing the programm.
  character(len=:), allocatable :: folder1, folder2
  character(len=50) :: arg1, arg2

  ! floating point arrays
  !   state     - state vector
  !   state_ref - reference state
  !   A_inv     - inverse Matrix needed for crank-nicolson-scheme
  real(dp), dimension(:), allocatable   :: state, state_ref
  real(dp), dimension(:,:), allocatable :: A_inv

  ! floating point values
  !   dt   - current time step
  real(dp) :: dt

  ! stat
  integer :: delstat

  ! dummy variable 
  integer :: i

  ! set folder
  if (iargc() == 1) then
    call getarg(1, arg1)
    folder1 = arg1
    folder2 = folder1

  else if (iargc() == 2) then
    call getarg(1, arg1)
    call getarg(2, arg2)
    folder1 = arg1
    folder2 = arg2

  else
    folder1 = fullname(".")
    folder2 = folder1
  end if

  !*******************
  ! CONVERGENCE TEST
  !*******************

  ! initialization
  write(*,*) '============================================================'
  write(*,*) ' Setting up Convergence Test:'
  write(*,*) '============================================================'
  !write(*,*) ''
  write(*,*) ' output folder (EOC):    ', color(trim(folder1),'yellow')
  write(*,*) ' output folder (states): ', color(trim(folder2),'yellow')

  call info( "  reading parameters from 'parameters.txt'")
  call read_param('parameters.txt', skip_first_char=.true.)
  call end_info()

  write(*,*) ''
  write(*,'(A,I13.1)')  '    initial amount of cells:     ', ncells
  write(*,'(A,I13.1)')  '    amount of ghost cells:       ', ngcells
  write(*,'(A,F13.4)')  '    maximal time:                ', tend
  write(*,'(A,F13.4)')  '    maximal timestep:            ', dt_max
  write(*,'(A,F13.4)')  '    cfl number:                  ', cfl
  write(*,'(A,F13.4)')  '    velocity:                    ', vel
  write(*,'(A,I13.1)')  '    reference resolution (cells) ', ncells_ref
  write(*,*)   '============================================================'

!  call info("  storing initial amount of cells")
!  ncells_ini = ncells
!  call end_info

!  ncells = ncells_ref
!  write(*,*)   ' setting up reference solution:'
!  write(*,*)   '============================================================'
!  write(*,*)   '     amount of cells:      ', ncells
!  write(*,*)   '============================================================'

!  call info('   adjusting time step')
!  dt_max = 1.0_dp/(3.0_dp)**(runs+1) * dt_max
!  call end_info()

!  ! initialize reference integration
!  call info( '   initializing state and setting parameters')
!  call initialize(state, state_ref, A_inv, dt)
!  call end_info()

!  write(*,*)   '  start integration:'
!  write(*,*)   ''

  ! integration
!  do while (time < tend)
!    call crank_nicolson_2(state, dt, A_inv, 1)
!    call show_progress()

!    time = time + dt
    
!    if ( write_cond(time, t_wo, 1.0d-16) ) then
!        call write_state(state, folder2)
!    end if
    
!  end do

!  call end_progress()

!  write(*,*) ''
!  write(*,*) color('   Done!','green')
!  write(*,*) ''

!  call info(  '  store reference frame')
!  state_ref = state
!  call end_info

!  call info(  '  reseting resolution and time step')
!  ncells = ncells_ini
!  dt_max = dt_max * 3.0_dp**(runs+1)
!  call end_info

  write(*,*) ''
  write(*,*) '  starting convergence test:'
  write(*,*) ''

  ! loop over resplutions
  do i = 1, runs

    ! improve resolution
    if (i > 1) then
      ncells = 		2 * ncells
      dt_max = 0.5_dp * dt_max
    end if

    write(*,*) '============================================================'
    write(*,*) '    amount of cells:       ', ncells
    write(*,*) '============================================================'

    if (allocated(state)) then
      call info( "  deallocating old arrays")
      deallocate(state, state_ref, A_inv)
      call end_info()
    end if


    call info( '  initializing state and set parameters')
    call initialize(state, state_ref, A_inv, dt)
    call end_info()

    write(*,*) ''
    write(*,*) ' starting integration:'
    
    ! integration
    call start_progress()

    do while (time < tend)
      call crank_nicolson_2(state, dt, A_inv, 1)
      call show_progress()

      time = time + dt
	 
      if ( write_cond(time, t_wo, 1.0d-16) ) then
        call write_state(state, folder2)
      end if
    end do

    if ( write_cond(time, t_wo, 1.0d-16) ) then
       call write_conv(state, state_ref, folder1)
    end if

    call end_progress()

    write(*,*) ''
    write(*,*) color('   Done!','green')
    write(*,*) ''

  end do

  call system('rm '//trim(folder1)//'*.bak', delstat)
  if (delstat /= 0) then
    call system('rm '//trim(folder1)//'/*.bak')
  end if


  write(*,*) '============================================================' 
  call info('  calculating order of convergence')
  call system('python src/add_eoc.py '//trim(folder1))
  call end_info()

  write(*,*) '============================================================'
  write(*,*) color('   CONVERGENCE TEST SUCCESSFULLY DONE!', 'green')
  write(*,*) '============================================================'

end program conv_test
