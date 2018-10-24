program test

  use progressbar
  implicit none

  integer :: i
  do i=1, 100
    call pbar(i, 100, 10)
  end do 
end program test

