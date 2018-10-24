module progressbar
  
  use colorize

  implicit none
  private

  public :: pbar

  contains

    ! methods

    subroutine pbar(it, max_it, nsteps, par_c, fill_c, per_c, t_c)
      ! Task:
      !  shows progress of the current task.
      ! Input:
      !  it     - current iteration,
      !  max_it - total number of iterations,
      !  nsteps - amount of steps to be shown in the bar,
      !  par_c  - [optional] color of the parenthesis,
      !  fill_c - [optional] color of the steps within the parenthesis,
      !  per_c  - [optional] color of the displayed percentages,
      !  t_c    - [optional] color of the displayed estimated remaining time.
      ! Output:
      !  progressbar str.
      implicit none
    
      integer,          intent(in)           :: it
      integer,          intent(in)           :: max_it
      integer,          intent(in)           :: nsteps
      character(len=*), intent(in), optional :: par_c 
      character(len=*), intent(in), optional :: fill_c
      character(len=*), intent(in), optional :: per_c
      character(len=*), intent(in), optional :: t_c
     
      integer :: per   ! actual progress in percentage 
      integer :: steps ! actual amount steps shown in the bar.
      integer :: hours, minutes, seconds 
      
      real :: tarray(2), current, remains

      integer :: i     ! dummy

      call etime(tarray, current)
      remains = current * (max_it / real(it) - 1.0) 
      hours = floor(remains / 3600)
      minutes = floor((remains - hours * 3600) / 60)
      seconds = nint(remains - (hours * 3600 + minutes * 60))

      per = nint( it/dble(max_it) * 100 ) 
      steps = floor(per * nsteps / 100.0)

      ! clear the line
      do i = 1, 20+2+nsteps+1+4+1+4+8
        write(6,'(a)',advance='no') '\b'
      end do

      if (per == 100) then 
        write(6, '(a)') ' -> Done.          ['
      else
        write(6,'(a)',advance='no') ' -> in progress... [' 
      end if


      if (steps .LE. 0) then
        do i = 1, nsteps
          write(6, '(a)', advance='no') ' '
        end do
      else if ( (steps .GT. 0) .and. (steps .LT. nsteps) ) then
        do i = 1, nsteps 
          write(6, '(a)', advance='no') '='
        end do
        write(6, '(a)', advance='no') '>'
      else 
        do i = 1, nsteps
          write(6, '(a)', advance='no') '='
        end do
      end if

      write(6, '(a)', advance='no') '] '
      write(6, '(I3.1)', advance='no') per
      write(6, '(a)', advance='no') '% ETA '
      write(6, '(I4.1)', advance='no') hours
      write(6, '(a)', advance='no') ':'
      write(6, '(I2.2)', advance='no') minutes
      write(6, '(a)', advance='no') ':' 
      write(6, '(I2.2)', advance='no') seconds

    end subroutine pbar
      

end module progressbar

