!******************************************************************************
! MODULE BASIC_IO: io for advection module
!******************************************************************************
module convergence_test_io
  use io_helpers
  use globals

  implicit none
  save

  ! integer values
  !   runs    - amount of runs: each run improves the resolution by a factor 2,
  !   scaling - scaling faktor for increasing resolution.
  integer :: runs
  integer :: scaling

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
          read(20,*) desc, diff_const
          read(20,*) desc, scaling
        ! skip header descripstion and read value line
        else
          read(20,*) ncells, ngcells, tend, dt_max, cfl, vel, runs, diff_const,&
            scaling
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


    subroutine write_conv(state, state_ref, folder, suffix_in)
      ! Task:
      !   writes out the number of cells 'ncells' and the error between initial
      !   state and actual state.
      ! Input:
      !   folder      -  output folder,
      !   state       -  current state,
      !   state_ref   -  reference state,
      !   suffix      -  optional, controls the data type [default = dat].

      implicit none

      character(len=:), allocatable, intent(in)           :: folder
      character(len=:), allocatable, intent(in), optional :: suffix_in
      character(len=4)                                    :: suffix

      ! dummy
      integer :: i

      ! state & reference state
      real(dp), dimension(ncells+2*ngcells), intent(in) :: state, state_ref

      ! error
      real(dp) :: err
	  real(dp) :: ratio1, ratio2

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
     
	    err = sum( abs(state(PHYmin:PHYmax) - state_ref(PHYmin:PHYmax)) )	
	    err = err / real(ncells, dp)

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

      write_cond = (mod(time,t_wo) < acc) .or. (1.0d0 - mod(time,t_wo) < acc)

    end function write_cond

end module convergence_test_io
