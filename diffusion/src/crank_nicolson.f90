!******************************************************************************
! MODULE CRANK_NICOLSON: implicit solver for hyperbolic PDEs
!******************************************************************************
module crank_nicolson
  use globals
  use basic_helpers, only: boundaries, update_dt

  implicit none
  private

  public :: crank_nicolson_2
  
  contains

    subroutine crank_nicolson_2(state, dt, A_inv, const_dt)
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
      real(dp), dimension(ncells, ncells), intent(in) :: A_inv

      ! floating point arrays
      !   B - matrix acting on actual state: B * u(t)
      !   state_inter - intermediate state:  B * u(t) = state_inter
      real(dp), dimension(ncells,ncells) :: B
      real(dp), dimension(ncells) :: state_inter

      ! dummy
      integer :: i

      ! lapack routine needed work arrays
      integer, dimension(ncells) :: ipiv
      integer :: info

      ! external procedures from lapack
      external DSYMV

      ! external DGESV  ! solves the system A*X=B for X, 
      !                 ! whereas A is symmertic n-by-n matrix
      !                 ! therefore using LU decomposition.
      !
      ! external DPOSV  ! like DGESV but A has to be hermitean and pos. def.
      !                 ! therefore using cholesky decomposition.

      ! update boundary conditions
      call boundaries(state)

      ! calculate new timestep
      call update_dt(dt, const_dt)

      ! update alpha   
      !alpha = diff_const * dt / (2.0_dp * dx**2)

      ! fill matrices
      !A(:,:) = 0.0_dp 
      !do i = 1, ncells-1
      !  A(i,i)   = 1.0_dp + 2.0_dp*alpha
      !  A(i,i+1) = -alpha
      !  A(i+1,i) = -alpha
      !end do

      !A(ncells, ncells) = 1.0_dp + 2.0_dp*alpha
      !A(1, ncells)      = -alpha
      !A(ncells, 1)      = -alpha

      B(:,:) = 0.0_dp
      do i = 1, ncells-1
        B(i,i) = 1.0_dp - 2.0_dp*alpha
        B(i,i+1) = alpha
        B(i+1,i) = alpha
      end do

      B(ncells, ncells) = 1.0_dp - 2.0_dp*alpha
      B(1, ncells) = alpha
      B(ncells, 1) = alpha
      

      state_inter(:) = 0.0_dp

      ! calculate B * u(t)
      !state_inter = matmul(B,state(PHYmin:PHYmax)) 
      call dsymv('U', size(state(PHYmin:PHYmax)), 1.0_dp, B, size(B,1), &
        state(PHYmin:PHYmax), 1, 0.0_dp, state_inter, 1)

      !call DGECON('I', size(A,1), A, ncells, sum(abs(A(1,:))), RCOND, WORK, IWORK, INFO) 

      !print*, rcond

      ! solve A * u(t+dt) = state_inter for u(t+dt)
      ! dgesv(  n   ,          - amount of linear equations
      !         nrhs,          - amount of colomns of inter_state
      !         A   ,          - matrix A
      !         LDA ,          - leading dimension of A
      !         ipiv,          - int array, for storing the pivot indices
      !         state_inter,   - intermediate state -> LATER SOLUTION
      !         LDstate_inter, - leading dimension of intermediate state
      !         info           - info integer
      !       )     
      !call DGESV(size(A,1), 1, A, size(A,1), ipiv, state_inter, size(A,1), info)
      !call DPOSV('U', ncells, 1, A, ncells, state_inter, ncells, info )

      ! check for the exact singularity.
      !if( info .gt. 0 ) then
      !   write(*,*)'The diagonal element of the triangular factor of A,'
      !   write(*,*)'U(',INFO,',',INFO,') is zero, so that'
      !   write(*,*)'A is singular; the solution could not be computed.'
      !   stop
      !end if
    
      ! update state
      !state(PHYmin:PHYmax) = matmul(A_inv,state_inter)

      call dsymv('U', size(state_inter), 1.0_dp, A_inv, size(A_inv,1), state_inter, &
        1, 0.0_dp, state(PHYmin:PHYmax), 1)

    end subroutine crank_nicolson_2
  
end module crank_nicolson
