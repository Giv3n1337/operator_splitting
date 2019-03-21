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
  !   ncells     - number of cells
  !   ngcells    - number of ghost cells
  !   PHYmin     - first index of physical domain
  !   PHYmax     - last index of physical domain
  !   ndims      - physical dimension
  integer :: ncells, ngcells, PHYmin, PHYmax, nDims

  ! floating point parameters
  !   pi - 3.1415......
  real(dp), parameter :: pi = 4.d0*datan(1.d0)

  ! floating point values
  !   time       - current time
  !   dt_max     - maximal time step
  !   tend       - maximal time
  !   t_wo       - write out time
  !   cfl        - cfl number
  !   dx         - cell size
  !   diff_const - diffusion constant
  real(dp) :: time, dt_max, tend, cfl, dx, t_wo, diff_const

  ! floating point arrays
  !   vel - fluid velocity
  !
  !   Coefficient Matrices: needed for crank nicolson procedure in order to 
  !                         solve u(t+dt) = A_inv * ( B * u(t))
  !     -> A_inv - inverse of coefficient matrix A
  !     -> B     - coefficient matrix B
  real(dp), dimension(:), allocatable :: vel  

  real(dp), dimension(:,:), allocatable :: A_inv, B

end module globals
