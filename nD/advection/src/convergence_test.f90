!******************************************************************************
! MAIN PROGRAM: convergence test of diffusion equation
!******************************************************************************
program convergence_test
  use convergence_test_io
  use initialization
  use runge_kutta

  implicit none

  ! characters
  !   folder - where to store output files,
  !   arg    - arguments passed in command line when executeing the programm.
  character(len=:), allocatable :: folder1, folder2
  character(len=50) :: arg1, arg2

  ! floating point arrays
  !   state     - state vector
  !   state_ref - reference state
  real(dp), dimension(:,:,:), allocatable :: state, state_ref

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
  write(*,*) '=============================================================='
  write(*,*) ' Setting up Convergence Test:'
  write(*,*) '=============================================================='
  !write(*,*) ''
  write(*,*) ' output folder (EOC):    ', color(trim(folder1),'yellow')
  write(*,*) ' output folder (states): ', color(trim(folder2),'yellow')

  call info( "  reading parameters from 'parameters.txt'")
  call read_param('parameters.txt', skip_first_char=.true.)
  call end_info()
  
  write(*,*) ''
  write(*,'(A,I1.1,A)')  '    Testing convergence of ', ndims,'D-Advection-Eq:'
  write(*,*) ''
  write(*,'(A,I13.1)')  '    initial amount of cells per dim:  ', ncells
  write(*,'(A,I13.1)')  '    amount of ghost cells:            ', ngcells
  write(*,'(A,F13.4)')  '    maximal time:                     ', tend
  write(*,'(A,F13.4)')  '    maximal timestep:                 ', dt_max
  write(*,'(A,F13.4)')  '    cfl number:                       ', cfl
  
  if (ndims == 2) then
    write(*,'(A,F13.4)')       '    velocity:                         ', vel(1)
    write(*,'(A,F13.4)')       '                                      ', vel(2)  
  else if (ndims == 3) then
    write(*,'(A,F13.4)')       '    velocity:                         ', vel(1)
    write(*,'(A,F13.4)')       '                                      ', vel(2) 
    write(*,'(A,F13.4)')       '                                      ', vel(3)
  end if
  
  write(*,'(A,F13.4)')  '    diffusion constant:               ', diff_const
  write(*,'(A,I13.1)')  '    max amount of cells:              ', ncells*scaling**(runs-1)
  write(*,*)   '=============================================================='


  write(*,*) ''
  write(*,*) '  starting convergence test:'
  write(*,*) ''

  ! loop over resolutions
  do i = 1, runs

    ! improve resolution
    if (i > 1) then
      ncells = scaling * ncells
      dt_max = 1.0_dp / real(scaling,dp) * dt_max
    end if

    write(*,*) '=============================================================='
    write(*,'(A,I15.2,A,I1.1)') '    amount of cells:              ', ncells,'^',ndims
    write(*,*) '=============================================================='

    if (allocated(state)) then
      call info( "  deallocating old arrays")
      deallocate(state, state_ref) !, A_inv, B)
      call end_info()
    end if

    call info( '  initializing state and set parameters')
    call initialize(state, state_ref, dt)
    !call update_coeff_matrices(dt)
    call end_info()

    !call info( '  write out initial state')
    !call write_state(state, folder2)
    !call end_info()

    write(*,*) ''
    write(*,*) ' starting integration:'
    
    ! integration
    call start_progress()

    do while (time < tend)
      !call strang_splitting(state, dt, 1)
      call rk2(state, dt, 1)
      call show_progress()

      time = time + dt
	 
      !if ( write_cond(time, 10.0_dp*dt_max, 1.0d-16) ) then
      !  call write_state(state, folder2)
      !end if
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


  write(*,*) '===================================================================' 
  call info('  calculating order of convergence')
  call system('python src/add_eoc.py '//trim(folder1))
  call end_info()

  write(*,*) '==================================================================='
  write(*,*) color('   CONVERGENCE TEST SUCCESSFULLY DONE!', 'green')
  write(*,*) '==================================================================='

end program convergence_test
