!******************************************************************************
! MODULE CRANK_NICOLSON: implicit solver for hyperbolic PDEs
!******************************************************************************
module crank_nicolson
  use globals
  use basic_helpers, only: boundaries, update_dt

  implicit none
  private

  public :: crank_nicolson_2
  public :: update_coeff_matrices

  ! TODO: write interface and SOR solver routine 

  contains

    subroutine crank_nicolson_2(state, dt, const_dt)
      ! Task: 
      !   2nd-order, implicit solver for hyperbolic PDEs -> u_t = D * u_xx
      !               => A * u(t+dt) = B * u(t)            | timestep t -> t+dt
      !               => u(t+dt) = A**(-1) * ( B * u(t) ) 
      !               => u(t+dt) = A**(-1) * u^(t)         | u^(t) = B * u(t)    
      ! Input:
      !   state    - actual state,
      !   dt       - previous time step,
      !   const_dt - indicator, whether dt = const. [0: false, 1: true]
      !   A_inv    - inverse matrix of A 
      ! Output:
      !   state - updated state vector,
      !   dt    - actual time step.

      implicit none

      integer, intent(in) :: const_dt

      real(dp), intent(inout) :: dt
      real(dp), dimension(ncells+2*ngcells), intent(inout) :: state

      ! floating point arrays
      !   state_inter - intermediate state:  B * u(t) = state_inter
      real(dp), dimension(ncells) :: state_inter

      ! dummy
      integer :: i
 
      ! lapack info
      integer :: info

      ! external procedures from lapack
      external DSYMV    ! matrix multiplication for symmetric matrices.
      
      ! update boundary conditions
      call boundaries(state)

      ! calculate new timestep
      call update_dt(dt, const_dt)

      if ( const_dt == 0 ) then
        call update_coeff_matrices(dt)
        ! call SOR(....)
        print*, 'not implemented jet!!!' !<----- TODO write SOR method
        return
      end if
      
      state_inter(:) = 0.0_dp

      ! calculate B * u(t)
      call dsymv('U', size(state(PHYmin:PHYmax)), 1.0_dp, B, size(B,1), &
        state(PHYmin:PHYmax), 1, 0.0_dp, state_inter, 1)

      ! calculate A^(-1)* (B * u(t)) = u(t+dt)
      call dsymv('U', size(state_inter), 1.0_dp, A_inv, size(A_inv,1), state_inter, &
        1, 0.0_dp, state(PHYmin:PHYmax), 1)

    end subroutine crank_nicolson_2
  
    
    subroutine update_coeff_matrices(dt)
      ! Task:
      !   updates the matrices A^(-1) and B depending on the time step.
      ! Input:
      !   dt - current time step.
      ! Output:
      !   A_inv - inverse of coefficent matrix A,
      !   B     - coefficent matrix B.  

      real(dp), intent(in) :: dt
      real(dp) :: alpha

      ! dummy
      integer :: i

      ! arrays needed by lapack routines
      real(dp), dimension(ncells, ncells) :: work
      integer, dimension(ncells) :: ipiv

      ! lapack info 
      integer :: info

      ! lapack routines
      external DSYTRI   ! symmetric matrix inversion
	    external DSYTRF   ! Bunch-Kaufman factorization of a symmetric matrix       

      ! update alpha   
      alpha = diff_const * dt / (2.0_dp * dx**2)

      ! fill matrix A and inverte it.
      A_inv(:,:) = 0.0_dp 
      do i = 1, ncells-1
        A_inv(i,i)   = 1.0_dp + 2.0_dp*alpha
        A_inv(i,i+1) = -alpha
        A_inv(i+1,i) = -alpha
      end do

      A_inv(ncells, ncells) = 1.0_dp + 2.0_dp*alpha
      A_inv(1, ncells)      = -alpha
      A_inv(ncells, 1)      = -alpha

      ! inversion
      call dsytrf('U', size(A_inv,1), A_inv, size(A_inv,1), ipiv, work, &
		    64*size(A_inv,1), info)
	  
	    if (info /= 0) then
		    stop 'matrix is numerically singular!'
	    end if

      call dsytri('U', size(A_inv,1), A_inv, size(A_inv,1), ipiv, work, &
		    info)

      ! fill matrix B
      B(:,:) = 0.0_dp
      do i = 1, ncells-1
        B(i,i) = 1.0_dp - 2.0_dp*alpha
        B(i,i+1) = alpha
        B(i+1,i) = alpha
      end do

      B(ncells, ncells) = 1.0_dp - 2.0_dp*alpha
      B(1, ncells) = alpha
      B(ncells, 1) = alpha

    end subroutine update_coeff_matrices

end module crank_nicolson
