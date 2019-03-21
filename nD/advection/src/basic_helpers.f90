!******************************************************************************
! MODULE BASIC_HELPERS: contains basic methods and interfaces
!******************************************************************************
module basic_helpers
  use globals

  implicit none
  private

  interface boundaries
    module procedure boundaries_2D
    module procedure boundaries_2D_nargs
    module procedure boundaries_3D_nargs
  end interface

  ! pulbic methods
  public :: update_dt     ! update_dt(dt, const_dt) 
  public :: boundaries    ! boundaries(array)

  contains

    subroutine boundaries_2D(state)
      ! Task:
      !   updates periodic boudary conditions of a 2-dim array.
      ! Input:
      !   state - state vector
      ! Output:
      !   state - updated state vector.

      implicit none

      real(dp), dimension(:,:), intent(inout) :: state

	  ! first dim
	  state(PHYmin:PHYmax,1:ngcells) = state(PHYmin:PHYmax,PHYmax+1-ngcells:PHYmax)
      state(PHYmin:PHYmax,PHYmax+1:) = state(PHYmin:PHYmax,PHYmin:PHYmin-1+ngcells)
	  
	  ! second dim
      state(1:ngcells,PHYmin:PHYmax) = state(PHYmax+1-ngcells:PHYmax,PHYmin:PHYmax)
      state(PHYmax+1:,PHYmin:PHymax) = state(PHYmin:PHYmin-1+ngcells,PHYmin:PHYmax)

    end subroutine boundaries_2D

    subroutine boundaries_2D_nargs(state)
      ! Task:
      !   updates periodic boudary conditions of a 2-dim array.
      ! Input:
      !   state - state vector
      ! Output:
      !   state - updated state vector.

      implicit none

      real(dp), dimension(:,:,:), intent(inout) :: state
	  integer :: i
	  
	  if (size(state,3) .eq. size(state,1)) then
		call boundaries_3D(state)
		return
	  end if
		
	  do i = 1, size(state,3)
		! first dim
		state(PHYmin:PHYmax,1:ngcells,i) = state(PHYmin:PHYmax,PHYmax+1-ngcells:PHYmax,i)
		state(PHYmin:PHYmax,PHYmax+1:,i) = state(PHYmin:PHYmax,PHYmin:PHYmin-1+ngcells,i)
	  
		! second dim
		state(1:ngcells,PHYmin:PHYmax,i) = state(PHYmax+1-ngcells:PHYmax,PHYmin:PHYmax,i)
		state(PHYmax+1:,PHYmin:PHymax,i) = state(PHYmin:PHYmin-1+ngcells,PHYmin:PHYmax,i)
	  end do

    end subroutine boundaries_2D_nargs


	subroutine boundaries_3D(state)
	  ! Task:
      !   updates periodic boudary conditions of a 3-dim array.
      ! Input:
      !   state - state vector
      ! Output:
      !   state - updated state vector.
	  
	  implicit none
	  
	  real(dp), dimension(:,:,:), intent(inout) :: state
		  
      ! first dim
	  state(PHYmin:PHYmax,PHYmin:PHYmax,1:ngcells) = &
		state(PHYmin:PHYmax,PHYmin:PHYmax,PHYmax+1-ngcells:PHYmax)
      
      state(PHYmin:PHYmax,PHYmin:PHYmax,PHYmax+1:) = &
		state(PHYmin:PHYmax,PHYmin:PHYmax,PHYmin:PHYmin-1+ngcells)
	  
	  ! second dim
      state(PHYmin:PHYmax,1:ngcells,PHYmin:PHYmax) = &
		state(PHYmin:PHYmax,PHYmax+1-ngcells:PHYmax,PHYmin:PHYmax)
      
      state(PHYmin:PHYmax,PHYmax+1:,PHYmin:PHymax) = &
		state(PHYmin:PHYmax,PHYmin:PHYmin-1+ngcells,PHYmin:PHYmax)

      ! third dim 
      state(1:ngcells,PHYmin:PHYmax,PHYmin:PHYmax) = &
		state(PHYmax+1-ngcells:PHYmax,PHYmin:PHYmax,PHYmin:PHYmax)
      
      state(PHYmax+1:,PHYmin:PHYmax,PHYmin:PHYmax) = &
		state(PHYmin:PHYmin-1+ngcells,PHYmin:PHYmax,PHYmin:PHYmax)
	  
	end subroutine boundaries_3D


    subroutine boundaries_3D_nargs(state)
	  ! Task:
      !   updates periodic boudary conditions of a 3-dim array.
      ! Input:
      !   state - state vector
      ! Output:
      !   state - updated state vector.
	  
	  implicit none
	  
	  real(dp), dimension(:,:,:,:), intent(inout) :: state
	  integer :: i	  
		  
	  if (size(state,4) > 9) then
	    stop
	  end if
		
	  do i = 1, size(state,4)
		  
		! first dim
		state(PHYmin:PHYmax,PHYmin:PHYmax,1:ngcells,i) = &
	      state(PHYmin:PHYmax,PHYmin:PHYmax,PHYmax+1-ngcells:PHYmax,i)
      
		state(PHYmin:PHYmax,PHYmin:PHYmax,PHYmax+1:,i) = &
		  state(PHYmin:PHYmax,PHYmin:PHYmax,PHYmin:PHYmin-1+ngcells,i)
	  
		! second dim
		state(PHYmin:PHYmax,1:ngcells,PHYmin:PHYmax,i) = &
		  state(PHYmin:PHYmax,PHYmax+1-ngcells:PHYmax,PHYmin:PHYmax,i)
      
		state(PHYmin:PHYmax,PHYmax+1:,PHYmin:PHymax,i) = &
		  state(PHYmin:PHYmax,PHYmin:PHYmin-1+ngcells,PHYmin:PHYmax,i)

		! third dim 
		state(1:ngcells,PHYmin:PHYmax,PHYmin:PHYmax,i) = &
		  state(PHYmax+1-ngcells:PHYmax,PHYmin:PHYmax,PHYmin:PHYmax,i)
      
        state(PHYmax+1:,PHYmin:PHYmax,PHYmin:PHYmax,i) = &
		  state(PHYmin:PHYmin-1+ngcells,PHYmin:PHYmax,PHYmin:PHYmax,i)
	  
	  end do
	end subroutine boundaries_3D_nargs


    subroutine update_dt(dt, const_dt)
      ! Task:
      !   generates new timestep depending on cfl criterion.
      !      -  if dt > dt_max       =>  dt = dt_max.
      !      -  if const_dt == 1     =>  dt = dt_max.
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
        dt = cfl * dx / maxval(vel)
      else 
        dt = dt_max
      end if

      if ( (time + dt > tend) .and. (time < tend) ) then
        dt = tend - time
      end if

    end subroutine update_dt
    

end module basic_helpers
