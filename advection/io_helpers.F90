module io_helpers
  
  use forbear, only: bar_object
  use FACE, only: color=>colorize, color_sample=>colors_samples
  implicit none

  ! exportet methods
  private
  public :: info
  public :: end_info
  public :: color
  public :: color_sample
  public :: bar_object
  public :: fullname

  contains

    subroutine info(str)
      ! Task:
      !  prints an info text.
      ! Input:
      !  str - info message.
      ! Output:
      !  str .
      implicit none
    
      character(len=*) , intent(in) :: str
    
      write(*,'(a)',advance='no') trim(str) // '... '   	

    end subroutine info


    subroutine end_info()
      ! Task:
      !  closes the info message, if the task is done.
      ! Output:
      !  print 'Done.' in green.
      implicit none
    
      print '(a)', color('Done.', color_fg='green') 

    end subroutine end_info
    
    character(len=255) function fullname(fname)
      ! Task:
      !  returns the full path to file/folder 'fname', if not already applied.
      ! Input:
      !  fname - string, that contains folder or .
      ! Output:
      !  fullname - string, that contains the full path name.
      implicit none
    

      character(len=*), intent(in) :: fname 
      character(len=128) :: prefix
      
      if (fname(1:1) .eq. '/') then
      	fullname = fname
      else if (fname(1:2) .eq. '~/') then
        call getenv('HOME', prefix)
        fullname = prefix(:lnblnk(prefix)) // fname(2:lnblnk(fname))
      else if (fname(1:1) .eq. '.') then
        call getcwd(prefix)
        fullname = prefix
      else
        call getcwd(prefix)
        fullname = prefix(:lnblnk(prefix)) // '/' // fname(:lnblnk(fname))
      end if

      return
    end function fullname


end module io_helpers
