!******************************************************************************
! MODULE RK2: runge kutta second order
!******************************************************************************
module runge_kutta
  use globals
  use basic_helpers, only: boundaries, update_dt
  use advection_rhs

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
