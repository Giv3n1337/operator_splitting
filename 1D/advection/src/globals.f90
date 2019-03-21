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
