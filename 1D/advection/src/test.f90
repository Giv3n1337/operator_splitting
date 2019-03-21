module m1
  use io_helpers  
  implicit none
  save

  type(bar_object) :: bar
  real(8) :: t

  contains

    subroutine initialize()
      use io_helpers
      implicit none 

      call bar%initialize(                  &
        bracket_left_string ='progress [',  &
        bracket_right_string=']',           &
        filled_char_string='=',             &
        empty_char_String=' ',              &
        add_progress_percent=.true.,        &
        max_value=dble(100.0))
    end subroutine initialize


    subroutine start()
    
      implicit none
    
      call bar%start()
    end subroutine start

    subroutine update()
    
      implicit none
    
      call bar%update(current=t)    
    end subroutine update   

    subroutine terminate()
    
      implicit none
    
      call bar%destroy
    
    end subroutine terminate
    

end module m1




program test
  use m1
  implicit none

  integer :: i,j
  real :: x
  t = 0.0d0
  
  call initialize() 
  call start()
  do i = 1, 100
    do j = 1,1000000
      x = exp(sqrt(real(j))*real(j))**2
    end do
    t = t + real(1,8)
    call update()
  end do

  call terminate()

end program test











