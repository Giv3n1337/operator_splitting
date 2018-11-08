
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
  !   time    - current time
  !   dt_max  - maximal time step
  !   tend    - maximal time
  !   t_wo    - write out time
  !   cfl     - cfl number
  !   vel     - fuild velocity
  !   dx      - cell size
  real(dp) :: time, dt_max, tend, cfl, vel, dx, t_wo

end module globals

!******************************************************************************
! MODULE ADVECTION_IO: io for advection module
!******************************************************************************
module advection_io
  use io_helpers
  use globals

  implicit none
  save

  ! integer values
  !   runs - amount of runs: each run improves the resolution by a factor 2,
  integer :: runs

  ! floating point values
  !   sig - standard deviation
  real(dp) :: sig 

  ! process bar
  type(bar_object) :: bar ! progress bar

  contains

    ! methods
    subroutine read_param(fname, hl, skip_first_char)
      ! Task:
      !   read initial values from file.
      ! Input:
      !   fname - file name,
      !   hl    - optional, amount of header lines to skip,
      !   skip_first_char - optional, logical parameter, if .true. skips first
      !                     char while reading a line. [default = .false.]
      ! Output:
      !  TODO: <-- write output param.

      implicit none

      character(len=*), intent(in)  :: fname
      integer, intent(in), optional :: hl
      logical, intent(in), optional :: skip_first_char

      ! io_stat
      integer :: open_stat, read_stat

      ! dummy character (description of a value to read in)
      character(len=8) :: desc

      ! dummy index
      integer :: i


      open(unit=20, file=fullname(fname), status='old', &
        action='read', form='formatted', iostat=open_stat)

      ! if file can be opened
      if (open_stat == 0) then

        ! skip heapper
        if (present(hl)) then
          do i=1, hl
            read(unit=20,fmt=*, iostat=read_stat)
            if (read_stat .ne. 0) then
              exit
            end if
          end do
        end if

        ! read in parameters TODO: assign right parameters
        if (present(skip_first_char) .and. skip_first_char .eqv. .true.) then
          read(20,*) desc, ncells
          read(20,*) desc, ngcells
          read(20,*) desc, tend
          read(20,*) desc, dt_max
          read(20,*) desc, cfl
          read(20,*) desc, vel
          read(20,*) desc, runs
          read(20,*) desc, sig
        ! skip header descripstion and read value line
        else
          read(20,*) ncells, ngcells, tend, dt_max, cfl, vel, runs, sig
        end if

      ! error opening file
      else
        write(*,*) 'ERROR: Cant open file. iostat: ', open_stat, '.'
      end if

      ! close file
      close(unit=20)

    end subroutine read_param



    subroutine write_state(state, folder, suffix_in)
      ! Task:
      !   writes actual state into folder 'folder'.
      ! Input:
      !   state   -  current state,
      !   folder  -  output folder,
      !   suffix  -  optional, controls data type [default = dat]

      implicit none

      real(dp), dimension(2*ngcells+ncells), intent(in) :: state

      character(len=:), allocatable, intent(in) :: folder
      character(len=4), optional,    intent(in) :: suffix_in


      ! io stat
      integer :: open_stat, write_stat
      logical :: ex

      ! filename
      character(len=255) :: fname
      character(len=7)   :: time_str
      character(len=4)   :: reso_str
      character(len=4)   :: suffix = 'txt '

      ! dummy
      integer :: i


      inquire(file = fullname(folder), exist=ex)
      if (.not. ex) then
        call system('mkdir -p '//fullname(folder))
      end if

      if (present(suffix_in)) then
        suffix = suffix_in
      end if

      write(time_str, '(f7.4)') time
      write(reso_str, '(I4)') ncells

      fname = trim(fullname(folder)) // '/' //'advection_n-' // &
        trim(adjustl(reso_str)) // '_t-' // trim(adjustl(time_str)) // &
        '.' // adjustl(suffix)

      ! open file
      open(unit=20, file=fname, status='replace', action='write', &
          form='formatted', iostat=open_stat)

      ! if file open
      if (open_stat == 0) then

        ! write out data
        do i = PHYmin, PHYmax
          write(unit=20,fmt=*, iostat=write_stat) state(i)

          ! if error or eof exit loop
          if (write_stat /= 0) then
            exit
          end if
        end do

        ! error writing data
        if (write_stat > 0) then
          write(*,*) 'ERROR: Unable to write file! iostat:', write_stat, '!'
        end if

      ! error opening the file
      else
        write(*,*) "ERROR: Can't open file! iostat: ", open_stat, '!'
      end if

      ! close file
      close(unit=20)


    end subroutine write_state


    subroutine write_conv(state, init_state, folder, use_steps, suffix_in)
      ! Task:
      !   writes out the number of cells 'ncells' and the error between initial
      !   state and actual state.
      ! Input:
      !   folder      -  output folder,
      !   state       -  current state,
      !   init_state  -  initital state,
      !   suffix      -  optional, controls the data type [default = dat].

      implicit none

      character(len=:), allocatable, intent(in)           :: folder
      character(len=:), allocatable, intent(in), optional :: suffix_in
      character(len=4)                                    :: suffix

      logical, intent(in), optional :: use_steps
      integer :: i

      real(dp), dimension(:), intent(in) :: init_state, state

      ! error
      real(dp) :: err

      ! io stat
      integer :: open_stat, write_stat
      logical :: ex

      ! filename
      character(len=:), allocatable :: fname
      character(len=7)   :: time_str

      suffix = 'txt '
      
      if (present(suffix_in)) then
        suffix = suffix_in
      end if

      inquire(file=fullname(folder), exist=ex)
      if (.not. ex) then
        call system('mkdir -p '//trim(fullname(folder)))
      end if

      write(time_str, '(f7.2)') time
      fname = trim(fullname(folder)) // 'advection_conv_test_' // 't-' &
        // trim(adjustl(time_str)) // '.' // adjustl(suffix)
      
      ! open file
      inquire(file=fname, exist=ex)
      if (ex) then
        open(unit=20, file=fname, status='old', position='append', &
          action='write', form='formatted', iostat=open_stat)
      else
        open(unit=20, file=fname, status='new', &
          action='write', form='formatted', iostat=open_stat)
      end if  
     
      err = sum(abs(state(PHYmin:PHYmax) - init_state(PHYmin:PHYmax))) / &
        real(ncells,dp)

      ! write resolution, error
      write(20,*,iostat=write_stat) ncells, err
 
      ! close file
      close(20)

      ! write backup file
      call backup(fname)

    end subroutine write_conv



    subroutine backup(fname)
      ! Task:
      !  creates a backup file incase the program crashes while writing out.
      ! Input:
      !  fname  -  filename.

      implicit none

      character(len=*), intent(in) :: fname
      character(len=:), allocatable :: copy_str

      logical :: ex

      inquire(file=fullname(fname), exist=ex)

      if (ex) then
        copy_str = 'cp ' // trim(fullname(fname)) // ' ' // &
          trim(fullname(fname)) // '.bak'
        call system(copy_str)
        return
      end if

      write(*,*) 'ERROR: file' // trim(fullname(fname)) // ' does not exist.'
      return
    end subroutine backup


    subroutine start_progress()
      ! Task:
      !   starts progress bar.
      implicit none

      call bar%start()
    end subroutine start_progress

    subroutine show_progress()
      ! Task:
      !   updates progress bar.
      implicit none

      call bar%update(current=time)
    end subroutine show_progress

    subroutine end_progress()
      ! Task:
      !   terminates the progress bar.
      implicit none
      call bar%update(current=time+0.0503_dp)
      call bar%destroy
    end subroutine end_progress


    logical function write_cond(time, t_wo, acc)
      ! Task: 
      !   checks whether write out condition is fullfilled.
      ! Input:
      !   time - current time,
      !   t_wo - time steps,
      !   acc  - accuracy.
      implicit none 

      real(dp), intent(in) :: time, t_wo, acc

      write_cond = mod(time,t_wo) < acc .or. 1.0d0 - mod(time,t_wo) < acc

    end function write_cond

end module advection_io



!******************************************************************************
! MODULE CONVERGENCE_HELPERS: io, initialization & dt routines
!******************************************************************************
module convergence_helpers
  use advection_io

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

    subroutine initialize(state, state_ini, dt)
      ! Task:
      !   allocates arrays, sets variables/parameters, set initial state
      ! Input:
      !   state     - state vector
      !   state_ini - initial state vector
      ! Outout:
      !   state     - state vector            -> fill
      !   state_ini - initial state vector    -> alloc
      !   dx        - grid size               -> calc
      !   dt        - time steps              -> set dt_max
      !   time      - current time            -> set 0.0
      !   t_wo      - write out time          -> calc
      !   PHYmax    - max index phys. domain  -> calc
      !   PHYmin    - min index phys. domain  -> calc
      
      implicit none

      ! floating point fields
      real(dp), dimension(:), allocatable, intent(inout) :: state, state_ini

      ! floating point values
      real(dp), intent(inout) :: dt

      ! dummy variable
      integer :: i

      ! calculate/set values/parameters and allocate arrays
      call initialize_param(state, state_ini, dt)

      ! set initial condition
      do i = 1, ncells
        state(ngcells+i) = 1.0_dp / (sig*sqrt(2.0_dp*pi)) * &     ! Gaussian
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

    subroutine initialize_param(state, state_ini, dt)
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
  use globals

  implicit none
  private

  ! pulbic methods
  public :: update_dt     ! update_dt(dt, const_dt) 
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



    subroutine update_dt(dt, const_dt)
      ! Task:
      !   generates new timestep depending on cfl criterion.
      !      -  if dt > dt_max       =>  dt = dt_max.
      !      -  if const_dt == 1     =>  dt = st_max.
      !      -  if time + dt > tend  =>  dt = tend - time.   
      ! Input:
      !   dt - last time step,
      !   const_dt - indicator: dynamic or constant time steps,
      !   last_dt  - logical to prevent endless loop.
      ! Output:
      !   dt - new time step.
    
      implicit none
    
        real(dp), intent(inout) :: dt
        integer,  intent(in)    :: const_dt
      
      if (const_dt == 0) then
        dt = cfl * dx / vel
      else 
        dt = dt_max
      end if

      if ( (time + dt > tend) .and. (time < tend) ) then
        dt = tend - time
      end if

    end subroutine update_dt
    

end module basic_helpers




!******************************************************************************
! MODULE RECONSTRUCTION: reconstruct interface values
!******************************************************************************
module reconstruction
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

end module reconstruction




!******************************************************************************
! MODULE ADVECTION_FLUX: flux calculation -> rhs
!******************************************************************************
module flux
  use globals
  use basic_helpers, only: boundaries
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

      real(dp), dimension(2*ngcells+ncells), intent(in)    :: state
      real(dp), dimension(2*ngcells+ncells), intent(inout) :: fluxes

      real(dp), dimension(2*ngcells+ncells,2) :: interfaces

      ! reconstruct interface values
      call recon_lin_van_leer(state, interfaces)
      !call godunov(state, interfaces) 

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
module right_hand_side
  use globals
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

      real(dp), dimension(ncells+2*ngcells) :: fluxes

      ! compute fluxes through interfaces
      call fluxes_adv_upwind(state, fluxes)

      ! calculating rhs
      res(:) = (fluxes(PHYmin:PHYmax) - fluxes(PHYmin+1:PHYmax+1)) / dx
    end subroutine rhs

end module right_hand_side




!******************************************************************************
! MODULE RK2: runge kutta second order
!******************************************************************************
module runge_kutta
  use globals
  use basic_helpers, only: boundaries, update_dt
  use right_hand_side

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
      call update_dt(dt, const_dt)
      call rhs(state, k2)
     

      k1(PHYmin:PHYmax) = state(PHYmin:PHYmax) + dt*k2

      call boundaries(k1)
      call rhs(k1, k2)

      state(PHYmin:PHYmax) = &
        0.5_dp * ( state(PHYmin:PHYmax) + k1(PHYmin:PHYmax) + dt*k2 )

    end subroutine rk2


end module runge_kutta





!******************************************************************************
! MAIN PROGRAM
!******************************************************************************
program conv_test
  use advection_io
  use convergence_helpers
  use runge_kutta

  implicit none

  ! characters
  !   folder - where to store output files,
  !   arg    - arguments passed in command line when executeing the programm.
  character(len=:), allocatable :: folder1, folder2
  character(len=50) :: arg1, arg2

  ! floating point arrays
  !   state     - state vector
  !   state_ini - initial state
  real(dp), dimension(:), allocatable :: state, state_ini

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

  call info( "  reading in parameters from 'parameters.txt'")
  call read_param('parameters.txt', skip_first_char=.true.)
  call end_info()

  write(*,'(A,I13.1)')  '    amount of cells:       ', ncells
  write(*,'(A,I13.1)')  '    amount of ghost cells: ', ngcells
  write(*,'(A,F13.4)')  '    maximal time:          ', tend
  write(*,'(A,F13.4)')  '    maximal timestep:      ', dt_max
  write(*,'(A,F13.4)')  '    cfl number:            ', cfl
  write(*,'(A,F13.4)')  '    velocity:              ', vel

  ! loop over resplutions
  do i = 1, runs

    ! improve resolution
    if (i > 1) then
      ncells =     2 * ncells
      dt_max = 0.5d0 * dt_max
    end if

    write(*,*) '============================================================'
    write(*,*) '    amount of cells:       ', ncells
    write(*,*) '============================================================'

    if (allocated(state)) then
      call info( "  deallocating old 'state' arrays")
      deallocate(state, state_ini)
      call end_info()
    end if

    call info( '  initializing state and set parameters')
    call initialize(state, state_ini, dt)
    call end_info()

    call info( '  copying initial state')
    state_ini(:) = state(:)
    call end_info()

    !call info('  write out initial state')
    !call write_state(state, folder2)
    !call end_info()
    
    write(*,*) ''
    write(*,*) ' start integration:'
    
    ! integration
    call start_progress()

    do while (time < tend)
      call rk2(state, dt, 1)
      call show_progress()

      time = time + dt

      
      if ( write_cond(time, dt_max, 1.0d-16) ) then
        call write_state(state, folder2)
      end if

      if ( write_cond(time, t_wo, 1.0d-16) ) then
         call write_conv(state, state_ini, folder1)
      end if
    end do

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
  write(*,*) color(' CONVERGENCE TEST SUCCESSFULLY DONE!', 'green')
  write(*,*) '============================================================'

end program conv_test
