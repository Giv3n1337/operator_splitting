!******************************************************************************
! MODULE BASIC_HELPERS: contains basic methods andinterfaces
!******************************************************************************
module basic_helpers
  use globals

  implicit none
  private

  ! pulbic methods
  public :: update_dt     ! update_dt(dt, const_dt) 
  public :: boundaries    ! boundaries(array)

  interface boundaries
    module procedure boundaries_1D
    module procedure boundaries_nD
  end interface

  contains

    subroutine boundaries_1D(state)
      ! Task:
      !   updates periodic boundary conditions of a 1D array.
      ! Input:
      !   state - state vector.
      ! Output:
      !   state - updated state vector.

      implicit none

      real(dp), dimension(:), intent(inout) :: state

      state(1:ngcells) = state(PHYmax+1-ngcells:PHYmax)
      state(PHYmax+1:) = state(PHYmin:PHYmin-1+ngcells)

    end subroutine boundaries_1D


    subroutine boundaries_nD(state)
      ! Task:
      !   updates periodic boudary conditions of a n-dim array.
      ! Input:
      !   state - state vector
      ! Output:
      !   state - updated state vector.

      implicit none

      real(dp), dimension(:,:), intent(inout) :: state
      integer :: i

      do i=1, size(state,2)
        state(1:ngcells,i) = state(PHYmax+1-ngcells:PHYmax,i)
        state(PHYmax+1:,i) = state(PHYmin:PHYmin-1+ngcells,i)
      end do

    end subroutine boundaries_nD



    subroutine update_dt(dt, const_dt)
      ! Task:
      !   generates new timestep depending on cfl criterion.
      !      -  if dt > dt_max       =>  dt = dt_max.
      !      -  if const_dt == 1     =>  dt = st_max.
      !      -  if time + dt > tend  =>  dt = tend - time.   
      ! Input:
      !   dt - last time step,
      !   const_dt - indicator: dynamic or constant time steps,
      !   last_dt  - logical to prevent endless loop.
      ! Output:
      !   dt - new time step.
    
      implicit none
    
        real(dp), intent(inout) :: dt
        integer,  intent(in)    :: const_dt
      
      if (const_dt == 0) then
        dt = cfl * dx / vel
      else 
        dt = dt_max
      end if

      if ( (time + dt > tend) .and. (time < tend) ) then
        dt = tend - time
      end if

    end subroutine update_dt
    

end module basic_helpers
