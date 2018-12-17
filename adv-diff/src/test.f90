program test

  implicit none

  ! solving the matrix equation A* x= b using LAPACK
  
  ! declarations, notice single precision
  double precision A(3,3), x(3), y(3), alpha, beta
  integer i, N, INCX, INCY, LDA

  double precision, dimension(3,3) :: work

  integer :: info

  integer, dimension(3) :: ipiv
 
  external DSYMV
  external DSYTRI   ! symmetric matrix inversion
	external DSYTRF   ! Bunch-Kaufman factorization of a symmetric matrix 
  
  ! define matrix A
  A(:,:) = 0.d0
  A(1,1)= 1.d0
  A(2,2)= 1.d0
  A(3,3)= 1.d0
  A(1,3)= 2.d0
  A(3,1)= 2.d0

  ! x
  X(1)= 3.d0
  X(2)= 4.d0
  X(3)= 5.d0
  
  Y(:) = 1.d0
  !Y(1) = 1.d0
  !Y(2) = 1.d0
  !Y(3) = 1.d0

  alpha = 1.0d0
  beta = 0.d0
  
  N = 3
  INCX = 1
  INCY = 1
  LDA = 3
  ! find the solution using the LAPACK routine SGESV
  !call DSYMV('U',N,alpha,A,LDA,X,INCX,beta,Y,INCY)
  call dsymv('U',size(X),1.d0,A,size(A,1),x,1,0.d0,y,1) 
  ! parameters in the order as they appear in the function call
  ! order of matrix A, number of right hand sides (b), matrix A,
  ! leading dimension of A, array that records pivoting,
  ! result vector b on entry, x on exit,leading dimension of b
  ! returnvalue   print the vector x
  do i= 1,3
      write(*,*) Y(i)
  end do

  Y = matmul(A,X)
  
  do i= 1,3
      write(*,*) Y(i)
  end do

  A(:,:) = 0.d0
  A(1,1) = 1.d0
  A(2,2) = 2.d0
  A(3,3) = 3.d0
  A(1,3) = 2.d0
  A(3,1) = 2.d0

  call dsytrf('U', size(A,1), A, size(A,1), ipiv, work, &
		64*size(A,1), info)

  call dsytri('U', size(A,1), A, size(A,1), ipiv, work, &
		info)
 
  do i= 1,3
    write(*,*) A(i,:)
  end do 

end program test
