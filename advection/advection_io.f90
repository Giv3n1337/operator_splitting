module advection_io
  use io_helpers
  use globals

  implicit none
  save
  
  type(bar_object) :: bar ! progress bar

  contains

    ! methods
    subroutine read_param(fname, hl, skip_first_char)
      ! Task:
      !  read initial values from file.
      ! Input:
      !  fname - file name,
      !  hl    - optional, amount of header lines to skip,
      !  skip_first_char - optional, logical parameter, if .true. skips first 
      !                    char while reading a line. [default = .false.] 
      ! Output:
      !  TODO: <-- write output param.
     
      implicit none
      
      character(len=*), intent(in)  :: fname
      integer, intent(in), optional :: hl
      logical, intent(in), optional :: skip_first_char

      ! io_stat
      integer :: open_stat, read_stat

      ! dummy character (description of a value to read in)
      character(len=:), allocatable :: desc  

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
        ! skip header descripstion and read value line  
        else 
          read(20,*) ncells, ngcells, tend, dt_max, cfl, vel
         end if
      
      ! error opening file
      else 
        write(*,*) 'ERROR: Cant open file. iostat: ', open_stat, '.'
      end if 
      
      ! close file 
      close(unit=20)  
      
    end subroutine read_param
    


    subroutine write_state(folder, suffix)
      ! Task:
      !  writes actual state into folder 'folder'.
      ! Input:
      !  folder  -  output folder,
      !  suffix  -  optional, controls data type [default = dat]
    	
      implicit none
    
      character(len=:), allocatable, intent(in)    :: folder
      character(len=4), optional,    intent(inout) :: suffix
      

      ! io stat     
      integer :: open_stat, write_stat
      logical :: ex

      ! filename
      character(len=255) :: fname
      character(len=7)   :: time_str
      character(len=4)   :: reso_str

      ! dummy
      integer :: i


      inquire(file = fullname(folder), exist=ex)
      if (.not. ex) then
        call system('mkdir '//fullname(folder))
      end if

      if (present(suffix)) then
        suffix = suffix
      else
        suffix = 'dat'
      end if

      write(time_str, '(f7.2)') time
      write(reso_str, '(I4)') ncells

      fname = fullname(folder) // 'advection_n=' // reso_str // '_t=' // & 
        time_str // '.' // trim(suffix)  

      ! open file 
      open(unit=20, file=fname, status='replace', action='write', &
          form='formatted', iostat=open_stat)
    
      ! if file open 
      if (open_stat == 0) then
        
        ! write out data
        do i = 1, ncells
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

    
    subroutine write_conv(folder, init_state, suffix)
      ! Task:
      !  writes out the number of cells 'ncells' and the error between initial
      !  state and actual state.
      ! Input:
      !  folder      -  output folder,
      !  init_state  -  initital state,  
      !  suffix      -  optional, controls the data type [default = dat].
      
      implicit none
    
      character(len=:), allocatable, intent(in)    :: folder
      character(len=4), optional,    intent(inout) :: suffix
      
      real(dp), dimension(ncells) :: init_state

      ! io stat     
      integer :: open_stat, write_stat
      logical :: ex

      ! filename
      character(len=255) :: fname
      character(len=7)   :: time_str

      ! dummy
      integer :: i

      if (present(suffix)) then
        suffix = trim(suffix)
      else
        suffix = 'dat'
      end if

      inquire(file = fullname(folder), exist=ex)
      if (.not. ex) then
        call system('mkdir '//fullname(folder))
      end if

      write(time_str, '(f7.2)') time

      fname = fullname(folder) // 'advection_conv_test' // '_t=' // time_str &
        // '.' // suffix  

      ! open file
      inquire(file=fname, exist=ex)
      if (ex) then 
        open(unit=20, file=fname, status='old', position='append', &
          action='write', form='formatted', iostat=open_stat)
      else 
        open(unit=20, file=fname, status='new', &
          action='write', form='formatted', iostat=open_stat)
      end if

      ! write resolution, error
      write(20,*,iostat=write_stat) ncells,  &
         sum(abs(state(1:ncells) - init_state)) / dble(ncells)
      
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
    

    subroutine show_progress()
      ! Task:
      !  updates progress bar.
      implicit none

      call bar%update(current=time)
    end subroutine show_progress

end module advection_io
