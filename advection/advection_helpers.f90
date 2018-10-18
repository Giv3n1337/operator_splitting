module advection_helpers
  use advection_io

  implicit none 
  save

  TODO: 


  subroutine cfl_crit(dt, a)                                                  
     ! Task:                                                                   
     !  calculate the next time step. (using global cfl number)                
     ! Input:                                                                  
     !  a  - array of current propagation speeds.                              
     ! Output:                                                                 
     !  dt - next time step.                                                   
                                                                               
     implicit none                                                             
                                                                               
     real(dp), intent(inout) :: dt                                             
     real(dp), dimension(:), intent(in) :: a                                   
                                                                               
     dt = cfl / sqrt( ((max(abs(a))) * ncells)**2 )                            
                                                                               
     if ( dt > dt_max ) then                                                   
       dt = dt_max                                                             
     end if                                                                    
                                                                               
     if ( time + dt > time_end ) then                                          
       dt = time_end - time                                                    
     end if                                                                    
  end subroutine cfl_crit

module advection_helpers                   
