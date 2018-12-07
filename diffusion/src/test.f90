program test

  implicit none

  ! solving the matrix equation A* x= b using LAPACK
  
  ! declarations, notice single precision
  DOUBLE PRECISION A(3,3), b(3)
  integer i, j, pivot(3), ok
  external DGESV
  
  ! define matrix A
  A(:,:) = 0.d0
  A(1,1)= 1
  A(2,2)= 2
  A(3,3)= 3
  
  ! rhs
  b(1)= 4.d0
  b(2)= 5.d0
  b(3)= 6.d0
  
  ! find the solution using the LAPACK routine SGESV
  call DGESV(3,1,A,3,pivot,b,3,ok)
  
  ! parameters in the order as they appear in the function call
  ! order of matrix A, number of right hand sides (b), matrix A,
  ! leading dimension of A, array that records pivoting,
  ! result vector b on entry, x on exit,leading dimension of b
  ! returnvalue   print the vector x
  do i= 1,3
      write(*,*) b(i)
  end do

end program test
