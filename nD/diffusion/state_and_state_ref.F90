     ! 10.d0 => arbitrary amplitude
     !  visc => diffusion konstant

     !decaying gauss  
     u(li,lj,lk,var) = 10.d0  * exp(-1.d0 * (&
                                ((i-0.5d0)*dx - 0.5d0)**2&
                              + ((j-0.5d0)*dy - 0.5d0)**2&
                              + ((k-0.5d0)*dz - 0.5d0)**2&
                                              ) / sig**2 &
                                      )

     u_ref(li,lj,lk,var) = 10.d0 / (1.d0 + 4.d0*tf*visc / sig**2)**(1.5d0)&
                              * exp(-1.d0 * (&
                                ((i-0.5d0)*dx - 0.5d0)**2&
                              + ((j-0.5d0)*dy - 0.5d0)**2&
                              + ((k-0.5d0)*dz - 0.5d0)**2&
                                              ) / (sig**2 + 4.d0*tf*visc)&
                                      )

     !decaying sin
     u(li,lj,lk,var) =  10.d0 * sin(2.d0*pi*((i-0.5d0)*dx&
                                            +(j-0.5d0)*dy&
                                            +(k-0.5d0)*dz))
     u_ref(li,lj,lk,var) =  10.d0 &
                              * exp(-1.d0*visc*(2.d0*pi)**2 *tf)   &
                              * sin(2.d0*pi*((i-0.5d0)*dx-vel(1)*tf&
                                            +(j-0.5d0)*dy-vel(2)*tf&
                                            +(k-0.5d0)*dz-vel(3)*tf))

