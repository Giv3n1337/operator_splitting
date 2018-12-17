!******************************************************************************
! MODULE LGS_SOLVER: contains routine for solving lgs: A*x = b,    (A = Matrix)
!******************************************************************************
module lgs_solver
  

  implicit none
  private

  public: SOR

    
  subroutine SOR(n,A,b,tollerance)
    ! Task: 
    !   iterative solver to estimate the solution of x for lgs like:
    !
    !        A*x = b, 
    !
    !   whereas x,b are vectors and A is a symmertric and pos.-def. Matrix, 
    !   which can be decomposed into:        A = D + U + L
    !     D - diagonal matrix,
    !     U - upper triangular matrix,
    !     L - lower triangular matri,
    ! Input:
    !   n - amount of equations
    !   A - coefficient matrix,
    !   b - right hand side,
    !   tollerance - difference between estimated solution of the last 
    !                iteration and the actual one.
    ! Output:
    !   b - solution of x.

    implicit none
   
    integer, intent(in) :: n

    real(dp), intent(in) :: tollerance
   
    real(dp), dimension(n,n), intent(in)  :: A
    real(dp), dimension(n), intent(inout) :: b

    ! work arrays
    real(dp), dimension(n,n) :: DL
    real(dp), dimension(n,n) :: U
    real(dp), dimension(n,n), intent :: x

    ! tuneable parameter for fast convergence
    real(dp) :: omega

    ! set initial solution as zero
    x(:,:) = 0.0_dp

    



  end subroutine SOR


end module lgs_solver
