!******************************************************************************
! MODULE: operator splitting
!******************************************************************************

module operator_splitting
  use globals
  use runge_kutta
  use crank_nicolson

  implicit none
  private

  public :: strang_splitting
  public :: update_coeff_matrices

  contains

    subroutine strang_splitting(state, dt, const_dt)
      ! Task:
      !  solving adv-diff eq using strang splitting (dt = const).
      ! Input:
      !  state - actual state vector,
      !  dt    - last time step,
      ! Output:
      !  state - updated state vector,
      !  dt    - actual time step.

      implicit none

      real(dp), dimension(ncells+2*ngcells), intent(inout) :: state

      real(dp), intent(inout) :: dt

      ! constant time step
      integer, intent(in) :: const_dt

      ! solve advection for t in [t, t + dt/2]
      call rk2(state, dt, const_dt)

      ! solve diffusion for t in [t, t + dt]
      call crank_nicolson_2(state, dt, const_dt)

      ! solve advection for [t + dt/2, t + dt]
      call rk2(state, dt, const_dt)

    end subroutine strang_splitting

end module
