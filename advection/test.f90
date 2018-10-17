program main
  use advection_helper
  implicit none

  character(len=:), allocatable :: full
  character(len=:), allocatable  :: system_str

  call read_param('parameters.txt', skip_first_char=.true., hl=1)
  print*, ncells, time_end, cfl

  call info(color('INFO: ', 'yellow') // 'this is a test info... ')
  call end_info()
  
  full = trim(fullname('parameters.txt'))
  system_str = 'cp ' // full // ' ' // full // '.bak'
  !call system(system_str)

  !print*, 'cp ' // trim(full) // ' ' // trim(full) // '.bak'
  !call system(system_str)
  ! print*, trim(system_str)
  call backup('parameters.txt')

end program main
