! ME634A end-sem assignment : Two-dimensional stretching of a circle in a shear velocity field using CLSVOF method.
! Code by Shreya Agrawal (160662) in Fortran 95

program main
        implicit none
        
        real :: Lx, Ly, del_t, tol, time, T, radius, center_x, center_y, dVf, resNR, tolNR, eps, min_d
        real :: F_low, F_up 
        real, dimension(:,:), allocatable :: u, v, F, Fstar, F_np1, phi, phi_star, &
                                                phi_np1, F_ul, F_ur, F_vb, F_vt, phi_aux, &
                                                phi_ul, phi_ur, phi_vb, phi_vt, nor_x, nor_y, d, s, d_prime 
        integer :: Nx, Ny, i, j, t_its, temp, initialF, itsNR, j_pp, i_pp, flag
        real, dimension(:), allocatable :: xe, xc, ye, yc, dxe, dxc, dye, dyc  ! grid edge and center coordinates for each grid-cell
        real, parameter :: pi = 3.1415927
        real, parameter :: one = 1.0


        ! Time-period for velocity
        T = 8

        ! Compuatational domain
        Lx = 1
        Ly = 1

        ! Cylinder radius and center
        radius   = 0.15
        center_x = 0.5
        center_y = 0.75
        
        ! Number of internal nodes in x,y directions
        Nx = 128
        Ny = 128

        ! Current time and iteration no. in time
        time  = 0
        t_its = 1

        min_d = 10000
        flag = 0

        F_low = 1.0e-10
        F_up = 1.0 - 1.0e-10
        
        temp = 0
        initialF = 0 ! Make =1 when initialising F (t = 0), and 0 at all other times 

        resNR = 1       ! Residual for NR iterations
        tolNR = 1.0e-5  ! tolerance for NR iterations
        itsNR = 0       ! NR iteration counter
        
        ! Thickness of the interface = epsilon
        eps = 1.75*(Lx/Nx)

        ! Matrix size allocation
        allocate(u(Nx+2, Ny+2))   ! x-velocity
        allocate(v(Nx+2, Ny+2))   ! y-velocity
        allocate(F(Nx+2, Ny+2))   ! VOF function
        allocate(phi(Nx+2, Ny+2))   ! LS function

        allocate(Fstar(Nx+2, Ny+2))   
        allocate(phi_star(Nx+2, Ny+2))

        allocate(F_np1(Nx+2, Ny+2))   ! VOF function at next time-step
        allocate(phi_np1(Nx+2, Ny+2)) ! LS function at next time-step

        allocate(phi_aux(Nx+2,Ny+2))

        allocate(F_ul(Nx+2, Ny+2))    ! VOF function at the cell edges
        allocate(F_ur(Nx+2, Ny+2))
        allocate(F_vb(Nx+2, Ny+2))
        allocate(F_vt(Nx+2, Ny+2))
        
        allocate(phi_ul(Nx+2, Ny+2))    ! LS function at the cell edges
        allocate(phi_ur(Nx+2, Ny+2))
        allocate(phi_vb(Nx+2, Ny+2))
        allocate(phi_vt(Nx+2, Ny+2))

        allocate(nor_x(Nx+2, Ny+2))
        allocate(nor_y(Nx+2, Ny+2))
   
        allocate(d(Nx+2, Ny+2))
        allocate(s(Nx+2, Ny+2))
        allocate(d_prime(Nx+2, Ny+2))

        allocate(xe(Nx+2))
        allocate(xc(Nx+2))
        allocate(ye(Ny+2))
        allocate(yc(Ny+2))
        
        allocate(dxe(Nx+2))
        allocate(dxc(Nx+2))
        allocate(dye(Ny+2))
        allocate(dyc(Ny+2))
        
        ! Initialisation
        u = 0
        v = 0
        F = 0
        phi = 0

        Fstar  = 0
        F_np1  = 0 

        phi_star  = 0
        phi_np1   = 0

        phi_aux = 10000

        F_ul = 0
        F_ur = 0
        F_vb = 0
        F_vt = 0
     
        phi_ul = 0
        phi_ur = 0
        phi_vb = 0
        phi_vt = 0
    
        nor_x = 0
        nor_y = 0
        
        d     = 0
        s     = 0
        d_prime = 0

        dxe = 0
        dxc = 0
        dye = 0
        dyc = 0
      
        ! Presently the grid is uniform, but can be made non-uniform through the following assignments
        do i = 1,Nx+2,1
                xe(i) = (i-1)*(Lx/Nx)
                xc(i) = (i-2+0.5)*(Lx/Nx)
        end do
        
        do j = 1,Ny+2,1
                ye(j) = (j-1)*(Ly/Ny)
                yc(j) = (j-2+0.5)*(Ly/Ny)
        end do

        do i = 2,Nx+2,1
               dxe(i) = xe(i) - xe(i-1)
        end do

        do i = 1,Nx+1,1
               dxc(i) = xc(i+1) - xc(i)
        end do

        do j = 2,Ny+2,1
               dye(j) = ye(j) - ye(j-1)
        end do

        do j = 1,Ny+1,1
               dyc(j) = yc(j+1) - yc(j)
        end do

        ! Ghost Cells grid-size
        dxe(1)    = dxe(2)
        dxc(Nx+2) = dxc(Nx+1)
        dye(1)    = dye(2)
        dyc(Ny+2) = dyc(Ny+1)
      
                ! velocity calculation for time = 0 
                
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                        u(i,j) = -( ( (SIN(pi*xe(i)))**2)*SIN(2*pi*yc(j))*COS((pi*time)/T) )
                                        v(i,j) =  ((SIN(pi*ye(j)))**2)*SIN(2*pi*xc(i))*COS((pi*time)/T)
                                end do
                        end do
                

        ! velocity BCs
        ! x = 0 (Inlet)
        u(1,:) = 0
        v(1,:) = -v(2,:)

        ! x = Lx (Outlet)
        u(Nx+1,:) = 0
        v(Nx+2,:) = -v(Nx+1,:)

        ! y = 0 (ground)
        v(:,1) = 0
        u(:,1) = -u(:,2)

        ! y = Ly (top)
        v(:,Ny+1) = 0
        u(:,Ny+2) = -u(:,Ny+1)

             
        ! LS function at t = 0 
        ! Pg.978 of Ref7
        ! Check this - phi<0 inside circle and phi>0 outside circle
        
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1
                phi(i,j) = -radius + SQRT((xc(i) - center_x)**2 + (yc(j) - center_y)**2)
              end do
           end do
        
        

        ! Boundary Conditions for LS function (Neumann)

        ! x = 0 (Inlet)
        phi(1,:) = phi(2,:)

        ! x = Lx (Outlet)
        phi(Nx+2,:) = phi(Nx+1,:)

        ! y = 0 (ground)
        phi(:,1) = phi(:,2)

        ! y = Ly (top)
        phi(:,Ny+2) = phi(:,Ny+1)

        
        ! Calculating VOF function at t = 0 
        initialF = 1
        call reconstruct(F, phi)

        ! Initialising F in all cells
        
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1

                if (phi(i,j) <= 0) then
                   F(i,j) = 0
                else
                   F(i,j) = 1
                end if

              end do
           end do
        

        ! Correcting F in interfacial cells
        
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1

                ! Finding cells where interface exists
                flag = 0
                
                  do j_pp = j-1,j+1,1
                    do i_pp = i-1,i+1,1
                       if ( ( (phi(i,j)*phi(i_pp,j_pp) < 0) &
                          .or. (phi(i,j) == 0) ) .and. ((i-i_pp)*(j-j_pp) <= 1) ) then

                   dVf = del_Vf(nor_x(i,j),nor_y(i,j),s(i,j),dxe(i),dye(j), phi(i,j))
                   F(i,j) = dVf/(dxe(i)*dye(j)) 

                        flag = 1

                      exit
                       else
                       end if

                    end do

                        if (flag == 1) then
                         exit
                        else
                        end if

                  end do
                
         

              end do
           end do
        

       

        ! Boundary Conditions for VOF function (Neumann)

        ! x = 0 (Inlet)
        F(1,:)   = F(2,:)

        ! x = Lx (Outlet)
        F(Nx+2,:)   = F(Nx+1,:)

        ! y = 0 (ground)
        F(:,1)   = F(:,2)

        ! y = Ly (top)
        F(:,Ny+2)   = F(:,Ny+1)

     
                ! Write to output files
                ! LS function
                open(unit = 100, file = 'phi_N128_t0.txt')
                do j = 1,Ny+2,1
                       
                        write(100,*) phi(:,j)
                end do
                close(100)

                ! VOF function
                open(unit = 200, file = 'F_N128_t0.txt')
                do j = 1,Ny+2,1
                        
                        write(200,*) F(:,j)
                end do
                close(200)

        

        initialF = 0 ! VOF initialisation is complete

        ! Solving in time
        do while ( time <= 8 )

                ! del_t calculation
                call delta_t(del_t)

                ! update present time
                time = time + del_t

                write(*,*) 'time =',time

                ! velocity calculation
                
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                        u(i,j) = -( ( (SIN(pi*xe(i)))**2 )*SIN(2*pi*yc(j))*COS((pi*time)/T) )
                                        v(i,j) =  ( (SIN(pi*ye(j)))**2)*SIN(2*pi*xc(i))*COS((pi*time)/T)                   
                                end do
                        end do
                              

        ! velocity BCs
        ! x = 0 (Inlet)
        u(1,:) = 0
        v(1,:) = -v(2,:)
       
        ! x = Lx (Outlet)
        u(Nx+1,:) = 0
        v(Nx+2,:) = -v(Nx+1,:)
        
        ! y = 0 (ground)
        v(:,1) = 0
        u(:,1) = -u(:,2)
        
        ! y = Ly (top)
        v(:,Ny+1) = 0
        u(:,Ny+2) = -u(:,Ny+1)

                     

        if ( MOD(t_its,2) == 0 ) then
        ! y-x sweep
                write(*,*) 'Going y-x!'

                call phiFaceCenter_ySweep(phi)
                call reconstruct(F,phi)
                call FfaceCenter_ySweep(F, phi)
                call calcStar(phi,phi_vb,phi_vt,phi_star,2)
                call calcStar(F,F_vb,F_vt,Fstar,2)

                call phiFaceCenter_xSweep(phi_star)
                call reconstruct(Fstar,phi_star)
                call FfaceCenter_xSweep(Fstar,phi_star)
                call calc_np1(phi_star,phi_ul,phi_ur,phi_np1,1)
                call calc_np1(Fstar,F_ul,F_ur,F_np1,1)


                  ! Truncating VOF Function
                   do j = 1,Ny+2,1
                      do i = 1,Nx+2,1

                        if ( (phi_np1(i,j) < -eps) .or. (F_np1(i,j) < F_low) ) then
                                F_np1(i,j) = 0
                        else if ( (phi_np1(i,j) > eps) .or. (F_np1(i,j) > F_up) ) then
                                F_np1(i,j) = 1
                        else
                        end if

                      end do
                   end do


                call reconstruct(F_np1,phi_np1)
                phi_aux = 10000
                call reinitialisePhi(F_np1,phi_np1)


        else
        ! x-y sweep
                write(*,*) 'Going x-y!'
                
                call phiFaceCenter_xSweep(phi)
                call reconstruct(F,phi)
                call FfaceCenter_xSweep(F,phi)
                call calcStar(phi,phi_ul,phi_ur,phi_star,1)
                call calcStar(F,F_ul,F_ur,Fstar,1)

                call phiFaceCenter_ySweep(phi_star)
                call reconstruct(Fstar,phi_star)
                call FfaceCenter_ySweep(Fstar,phi_star)
                call calc_np1(phi_star,phi_vb,phi_vt,phi_np1,2)
                call calc_np1(Fstar,F_vb,F_vt,F_np1,2)

         
                  ! Truncating VOF Function
                   do j = 1,Ny+2,1
                      do i = 1,Nx+2,1

                        if ( (phi_np1(i,j) < -eps) .or. (F_np1(i,j) < F_low) ) then
                                F_np1(i,j) = 0
                        else if ( (phi_np1(i,j) > eps) .or. (F_np1(i,j) > F_up) ) then
                                F_np1(i,j) = 1
                        else
                        end if

                      end do
                   end do


                call reconstruct(F_np1,phi_np1)
                phi_aux = 10000
                call reinitialisePhi(F_np1,phi_np1)

        end if
  
                phi = phi_np1
                F   = F_np1

                t_its = t_its + 1

                ! Write to output files

           if ( (time > 1.99) .and. (time < 2.01) ) then

                ! LS function
                open(unit = 300, file = 'phi_N128_t2.txt')
                do j = 1,Ny+2,1
                        
                        write(300,*) phi(:,j)
                end do
                close(300)

                ! VOF function
                open(unit = 400, file = 'F_N128_t2.txt')
                do j = 1,Ny+2,1
                        
                        write(400,*) F(:,j)
                end do
                close(400)
             
           else if ( (time > 3.94) .and. (time < 4.05) ) then

                ! Write to output files
                ! LS function
                open(unit = 500, file = 'phi_N128_t4.txt')
                do j = 1,Ny+2,1
                        
                        write(500,*) phi(:,j)
                end do
                close(500)

                ! VOF function
                open(unit = 600, file = 'F_N128_t4.txt')
                do j = 1,Ny+2,1
                       
                        write(600,*) F(:,j)
                end do
                close(600)

           else if ( (time > 5.99) .and. (time < 6.01) ) then

                ! Write to output files
                ! LS function
                open(unit = 700, file = 'phi_N128_t6.txt')
                do j = 1,Ny+2,1
                        
                        write(700,*) phi(:,j)
                end do
                close(700)

                ! VOF function
                open(unit = 800, file = 'F_N128_t6.txt')
                do j = 1,Ny+2,1
                        
                        write(800,*) F(:,j)
                end do
                close(800)

           else
           end if


        end do


                ! Write to output files
                ! LS function
                open(unit = 900, file = 'phi_N128_t8.txt')
                do j = 1,Ny+2,1
                        
                        write(900,*) phi(:,j)
                end do
                close(900)

                ! VOF function
                open(unit = 1000, file = 'F_N128_t8.txt')
                do j = 1,Ny+2,1
                        
                        write(1000,*) F(:,j)
                end do
                close(1000)


contains

        ! To reconstruct the interface using phi and F
        ! Gives nor_x, nor_y, nor_z, d (without mass conservation) and s (ensuring mass conservation)
        ! F_temp = F/Fstar/F_np1, phi_temp = phi/phi_star/phi_np1
        subroutine reconstruct(F_temp, phi_temp)

        implicit none
        integer :: i, j, i_prime, j_prime, i_p, j_p
        real :: tt, Fc, dx1, dy1, t1, t2, sc, s1, sm, w8, Hp, Xt, Yt 
        real, dimension(Nx+2, Ny+2), intent(in) :: F_temp, phi_temp
        real, dimension(3,3) :: A
        real, dimension(3)   :: x, b        

        tt = 0 ! use tt to normalise the plane 
        
        
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1
                            
                 ! Assign zero normals to all cells and update for only
                 ! interfacial cells
                 nor_x(i,j) = 0
                 nor_y(i,j) = 0
                 d(i,j)     = 0
                 s(i,j)     = 0

                 flag = 0

                 ! Finding cells where interface exists        
                 if ( ((F_temp(i,j) > F_low) .and. (F_temp(i,j) < F_up)) .or. (initialF == 1)  ) then
                    
                      do j_prime = j-1,j+1,1
                        do i_prime = i-1,i+1,1
                           if ( ( (phi_temp(i,j)*phi_temp(i_prime,j_prime) < 0) &
                             .or. (phi_temp(i,j) == 0) ) & 
                                .and. ((i-i_prime)*(j-j_prime) <= 1) ) then


                                   flag = 1
                                                                                                                
                                   ! Reconstruction using phi 

                                   A = 0
                                   b = 0
                                   x = 0
        
                                   
                                      do j_p = j-1,j+1,1
                                         do i_p = i-1,i+1,1

                                            Xt = xc(i_p) - xc(i)
                                            Yt = yc(j_p) - yc(j)

                                            if ((i_p == i) .and. (j_p == j)) then
                                                w8 = 52
                                            else
                                                w8 = 1
                                            end if
                                                
                                            if ( ABS(phi_temp(i_p,j_p)) <= eps ) then
                                                Hp = (1/(2*eps))*(1+COS((pi*phi_temp(i_p,j_p))/eps))
                                            else
                                                Hp = 0
                                            end if    
                                       

                                            A(1,1) = A(1,1) + w8*Hp*Xt*Xt
                                            A(2,1) = A(2,1) + w8*Hp*Xt*Yt
                                            A(3,1) = A(3,1) + w8*Hp*Xt

                                            A(1,2) = A(1,2) + w8*Hp*Xt*Yt
                                            A(2,2) = A(2,2) + w8*Hp*Yt*Yt
                                            A(3,2) = A(3,2) + w8*Hp*Yt

                                            A(1,3) = A(1,3) + w8*Hp*Xt
                                            A(2,3) = A(2,3) + w8*Hp*Yt
                                            A(3,3) = A(3,3) + w8*Hp

                                            b(1) = b(1) + w8*Hp*phi_temp(i_p,j_p)*Xt
                                            b(2) = b(2) + w8*Hp*phi_temp(i_p,j_p)*Yt
                                            b(3) = b(3) + w8*Hp*phi_temp(i_p,j_p)
                                                
                                         end do
                                      end do
                                                                                                                          
                                   ! Use Gauss Elimination for solving the matrix 
                                   call GaussElimination(A,x,b,3)

                                   ! Normalise x
                                   tt = ( x(1)**2 + x(2)**2 )**0.5
                                   nor_x(i,j) = x(1)/tt
                                   nor_y(i,j) = x(2)/tt
                                   d(i,j)     = x(3)/tt
                                   

                                   ! Calculate s using d - Pg.11 of CLSVOF doc       
                                   s(i,j) = d(i,j) + 0.5*( ABS(nor_x(i,j))*dxe(i) + &
                                              ABS(nor_y(i,j))*dye(j) )     

                                   ! Correcting s using F
                                   if (initialF == 0) then
                                         
                                       t1 = dxe(i)*ABS(nor_x(i,j))
                                       t2 = dye(j)*ABS(nor_y(i,j))

                                       ! Ensuring that dx1t >= dy1t
                                       dx1 = max(t1,t2)
                                       dy1 = min(t1,t2)

                                       sm = dx1 + dy1 
                                       sc = min(s(i,j),sm-s(i,j)) 
                                       Fc = min(F_temp(i,j),1-F_temp(i,j)) 
                                        
                                       sc = (2*Fc*dx1*dy1)**0.5
                                       s1 = sc
 
                                       if (sc <= dy1) then
                                         if (F_temp(i,j) <= 0.5) then
                                             s(i,j) = sc
                                         else 
                                             s(i,j) = sm - sc
                                         end if
                                       
					else 
                                         
                                         sc = Fc*dx1 + 0.5*dy1    
                                                
                                         if (F_temp(i,j) <= 0.5) then
                                              s(i,j) = sc
                                         else
                                              s(i,j) = sm - sc
                                         end if

                                       end if
                        
                                   else
                                   end if           
                                            
                                   exit

                           else
                           end if
                        end do

                         if (flag == 1) then
                         exit
                        else
                        end if


                      end do
                    
                                     
                 else
                 end if

              end do
           end do
                

        end subroutine reconstruct



        ! Returns the volume of fluid in a given region (sub_x, sub_y) within a cell
        ! The region must contain origin of s 
        function del_Vf(nor_xt, nor_yt, st, sub_x, sub_y,phi_ij)

        implicit none
        real :: nor_xt, nor_yt, st, sub_x, sub_y, del_Vf
        real :: t1, t2, sm, sc, Fc, dx1t, dy1t, phi_ij

        t1 = sub_x*ABS(nor_xt)
        t2 = sub_y*ABS(nor_yt)

        ! Ensuring that dx1t >= dy1t        
        dx1t = max(t1,t2)
        dy1t = min(t1,t2)


        sm = dx1t + dy1t
        sc = min(st,sm-st) ! st = s(i,j)

        ! For non-interfacial cells
        if (((nor_xt == 0.0) .and. (nor_yt == 0.0)) ) then

                if (phi_ij > 0) then
                        Fc = 1
                        del_Vf = Fc*sub_x*sub_y
                else
                        Fc = 0
                        del_Vf = Fc*sub_x*sub_y
                end if

        else
        ! For interfacial cells

        if ( sm <= st ) then
                Fc = 1
                del_Vf = Fc*sub_x*sub_y
        else if ( st < 0 ) then
                Fc = 0
                del_Vf = Fc*sub_x*sub_y
        else 
                if ( sc <= dy1t ) then

                        Fc = (sc**2)/(2*dx1t*dy1t) 

                else 

                        Fc = (sc**2 - (sc - dy1t)**2)/(2*dx1t*dy1t)

                end if

        
                if (st <= 0.5*sm) then
                      del_Vf = Fc*sub_x*sub_y 
                else
                      del_Vf = (1-Fc)*sub_x*sub_y
                end if

        end if

        end if

        end function del_Vf
        


        ! To calculate F on i-1/2 and i+1/2 faces
        ! F_sweep = F (x-y sweep)/Fstar (y-x sweep)
        subroutine FfaceCenter_xSweep(F_sweep,phi_sweep)

        implicit none
        integer :: i,j,k
        real, dimension(Nx+2, Ny+2),intent(in) :: F_sweep, phi_sweep
        real :: delta_Vf, ut, nxt, nyt, s_t
        ! ut = u(i,j), nxt = nor_x(i,j), nyt = nor_y(i,j), s_t = s(i,j)        

        
           do j = 2,Ny+1,1
              do i = 1,Nx+1,1

                ut = u(i,j)

                ! right face
                if (ut > 0) then
        
                   nxt = nor_x(i,j)
                   nyt = nor_y(i,j)
                   s_t = s(i,j)
                
                   if ( nxt*ut >= 0) then
                     delta_Vf = del_Vf(nxt,nyt,s_t,del_t*ABS(ut),dye(j), phi_sweep(i,j))                                        
                   else 
                     delta_Vf = F_sweep(i,j)*dxe(i)*dye(j) - &
                                del_Vf(nxt,nyt,s_t,dxe(i) - del_t*ABS(ut),dye(j),phi_sweep(i,j))
                   end if                
                        
                else if (ut < 0) then

                   nxt = nor_x(i+1,j)
                   nyt = nor_y(i+1,j)
                   s_t = s(i+1,j)

                   if ( nxt*ut >= 0) then
                     delta_Vf = del_Vf(nxt,nyt,s_t,del_t*ABS(ut),dye(j),phi_sweep(i+1,j))                                  
                   else
                     delta_Vf = F_sweep(i+1,j)*dxe(i+1)*dye(j) - &
                                del_Vf(nxt,nyt,s_t,dxe(i+1) - del_t*ABS(ut),dye(j),phi_sweep(i+1,j))   
                   end if
                                
                else                
                end if                                

                if (ut == 0) then
                        if (i == 1) then
                          F_ur(i,j) = F_sweep(i,j) + (dxe(i)/2)*((F_sweep(i+1,j) - &
                                          F_sweep(i,j))/(xc(i+1) - xc(i)))
                        else
                          F_ur(i,j) = F_sweep(i,j) + (dxe(i)/2)*((F_sweep(i+1,j) - &
                                          F_sweep(i-1,j))/(xc(i+1) - xc(i-1)))
                        end if

                else
                        F_ur(i,j) = delta_Vf/(ABS(ut)*del_t*dye(j))
                end if

                ! left face for next cell
                F_ul(i+1,j) = F_ur(i,j)

              end do
           end do
        

        end subroutine FfaceCenter_xSweep


        ! To calculate F on j-1/2 and j+1/2 faces
        ! F_sweep = F (y-x sweep)/Fstar (x-y zweep)
        subroutine FfaceCenter_ySweep(F_sweep,phi_sweep)

        implicit none
        integer :: i,j
        real, dimension(Nx+2, Ny+2),intent(in) :: F_sweep,phi_sweep
        real :: delta_Vf, vt, nxt, nyt, s_t
        ! vt = v(i,j), nxt = nor_x(i,j), nyt = nor_y(i,j), s_t = s(i,j)

        
           do j = 1,Ny+1,1
              do i = 2,Nx+1,1

                vt = v(i,j)

                ! top face
                if (vt > 0) then

                   nxt = nor_x(i,j)
                   nyt = nor_y(i,j)
                   s_t = s(i,j)

                   if ( nyt*vt >= 0) then
                     delta_Vf = del_Vf(nxt,nyt,s_t,dxe(i),del_t*ABS(vt),phi_sweep(i,j))
                   else
                     delta_Vf = F_sweep(i,j)*dxe(i)*dye(j) - &
                                del_Vf(nxt,nyt,s_t,dxe(i),dye(j) - del_t*ABS(vt),phi_sweep(i,j))
                   end if

                else if (vt < 0) then

                   nxt = nor_x(i,j+1)
                   nyt = nor_y(i,j+1)
                   s_t = s(i,j+1)

                   if ( nyt*vt >= 0) then
                     delta_Vf = del_Vf(nxt,nyt,s_t,dxe(i),del_t*ABS(vt),phi_sweep(i,j+1))
                   else
                     delta_Vf = F_sweep(i,j+1)*dxe(i)*dye(j+1) - &
                                del_Vf(nxt,nyt,s_t,dxe(i),dye(j+1) - del_t*ABS(vt),phi_sweep(i,j+1))
                   end if

                else
                end if

                if (vt == 0) then
                        if (j == 1) then
                          F_vt(i,j) = F_sweep(i,j) + (dye(j)/2)*((F_sweep(i,j+1) - &
                                          F_sweep(i,j))/(yc(j+1) - yc(j)))
                        else
                          F_vt(i,j) = F_sweep(i,j) + (dye(j)/2)*((F_sweep(i,j+1) - &
                                          F_sweep(i,j-1))/(yc(j+1) - yc(j-1)))
                        end if

                else
                        F_vt(i,j) = delta_Vf/(ABS(vt)*del_t*dxe(i))
                end if

                ! bottom face for the next cell
                F_vb(i,j+1) = F_vt(i,j)

              end do
           end do
       

        end subroutine FfaceCenter_ySweep


        ! Reinitialisation of LS function 
        ! F_tem = F_np1, phi_tem = phi_np1
        ! Corrects phi_tem according to the reconstructed interface
        subroutine reinitialisePhi(F_tem, phi_tem)

        real, dimension(Nx+2, Ny+2) :: d_cor, phi_tt
        real, dimension(Nx+2, Ny+2), intent(in) :: F_tem
        real, dimension(Nx+2, Ny+2), intent(out) :: phi_tem
        integer :: i,j,Kk,i_prime,j_prime, i_p,j_p,l,m
        integer :: inew, jnew, iin, jjn
        real, dimension(2) :: xv,xp,x_prime, ps1, ps2, pis, pr1, pr2, tmp1,tmp2
        real :: d_prod, nnx, nny, nnz, d_max, nxt, nyt, nzt, nxc, nyc, nzc, dt, dc, nor, var1, var2
       
        phi_aux  = 10000 
        phi_tt = 0
        min_d = 10000
        ! phi is determined accurately in (i-Kk,j-Kk) to (i+Kk,j+Kk) region surrounding interfacial cell (i,j,k)
        Kk = 4
        d_cor = 0
        d_prime = 0


        tmp1(1) = xe(1)
        tmp1(2) = ye(1)
        tmp2(1) = xe(2)
        tmp2(2) = ye(2)

        ! Maximum diagonal length of all computational cells
        d_max = Euc(tmp1,tmp2)
        
        ! d_cor is the corrected d obtained from the corrected s
        
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1
                d_cor(i,j) = s(i,j) - 0.5*( ABS(nor_x(i,j))*dxe(i) + &
                               ABS(nor_y(i,j))*dye(j) )
              end do
           end do
        

        ! d_prime is the distance between the reconstructed interface and origin (0,0)
        
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1
                d_prime(i,j) = d_cor(i,j) - ( nor_x(i,j)*xc(i) + &
                                 nor_y(i,j)*yc(j) )
              end do
           end do
       

        
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1

                flag = 0

                 ! Finding cells where interface exists
                 if ( (F_tem(i,j) > F_low) .and. (F_tem(i,j) < F_up) ) then
                    
                      do j_prime = j-1,j+1,1
                        do i_prime = i-1,i+1,1
                           if ( ( (phi_tem(i,j)*phi_tem(i_prime,j_prime) < 0) &
                              .or. (phi_tem(i,j) == 0) ) .and. ((i-i_prime)*(j-j_prime) <= 1) ) then

                                flag = 1
                        
                                nnx = nor_x(i,j)
                                nny = nor_y(i,j)

                                
                                  do j_p = j-Kk,j+Kk,1
                                    do i_p = i-Kk,i+Kk,1

                                  ! Checking that the cells near the interfacial cell (within a band of Kk) indeed are inside the
                                  ! computational domain  
                                  if ( (i_p >= 1) .and. (i_p <= Nx+2) .and. (j_p >= 1) .and. &
                                       (j_p <= Ny+2) ) then  
                                        
                                        ! x_prime = cell center (i_p,j_p) vector
                                        x_prime(1) = xc(i_p)
                                        x_prime(2) = yc(j_p)

                                        ! when i_p = i, j_p = j(interfacial cell)
                                        if ((i_p == i) .and. (j_p == j)) then
                                           phi_aux(i_p,j_p) = SIGN(one,F_tem(i_p,j_p)-0.5)*ABS(Dis(i,j,x_prime))
                                        else
                                          
                                           ! Determining xv = point on the  boundary of cell (i,j) with the shortest distance to cell center (i_p,j_p)
                                           l = max(-1, min(1, i_p - i))
                                           m = max(-1, min(1, j_p - j))

                                           if (l == -1) then
                                                xv(1) = xe(i-1)
                                           else if (l == 0) then
                                                xv(1) = xc(i)
                                           else if (l == 1) then
                                                xv(1) = xe(i)
                                           else
                                           end if

                                           if (m == -1) then
                                                xv(2) = ye(j-1)
                                           else if (m == 0) then
                                                xv(2) = yc(j)
                                           else if (m == 1) then
                                                xv(2) = ye(j)
                                           else
                                           end if

                                           if ( Dis(i,j,xv)*SIGN(one,F_tem(i_p,j_p)-0.5) <= 0 ) then
                                              phi_aux(i_p,j_p) = &
                                              SIGN(one,F_tem(i_p,j_p)-0.5)*min( Euc(xv,x_prime),ABS(phi_aux(i_p,j_p)) )
                                           else
                                           end if
                                                
                                              ! xp = projection of x_prime onto the interface in cell (i,j,k)
                                              xp(1) = x_prime(1) - Dis(i,j,x_prime)*nnx  
                                              xp(2) = x_prime(2) - Dis(i,j,x_prime)*nny

                                              ! Checking if xp is inside cell (i,j,k)
                                              if (((xp(1) > xe(i-1)) .and. (xp(1) < xe(i))) .and. &
                                                 ( (xp(2) > ye(j-1)) .and. (xp(2) < ye(j))) ) then

                                                phi_aux(i_p,j_p) = &
                                                SIGN(one,F_tem(i_p,j_p)-0.5)*min( Euc(xp,x_prime),ABS(phi_aux(i_p,j_p)) )

                                              else
                                    
                                                ! (Step 4)
                                                min_d = 10000
                                                
                                                 do jnew = j-1,j,1
                                                  do inew = i-1,i,1

                                                    
                                                     do jjn = j-1,j,1
                                                      do iin = i-1,i,1

                                                        ps1(1) = xe(inew)
                                                        ps1(2) = ye(jnew)
                                                        ps2(1) = xe(iin)
                                                        ps2(2) = ye(jjn)

                                                        d_prod = Dis(i,j,ps1)*Dis(i,j,ps2)
                                                            

                                                        ! edge corners separated by dye(j)         
                                                        if ( (jnew < jjn) .and. (inew == iin)) then
                                                          if (d_prod > 0) then
                                                          else
                                                            pis(1) = ps1(1)    
                                                            pis(2) = (nnx*pis(1) + d_prime(i,j))/(-nny)

                                                            min_d = min( Euc(pis,x_prime),ABS(min_d) )

                                                          end if
                                                         
                                                        ! edge corners separated by dxe(i)
                                                        else if ( (jnew == jjn) .and. (inew < iin)) then
                                                          if (d_prod > 0) then
                                                          else
                                                            pis(2) = ps1(2)
                                                            pis(1) = (nny*pis(2) + d_prime(i,j))/(-nnx)

                                                            min_d = min( Euc(pis,x_prime),ABS(min_d) )
       
                                                          end if
                                                        else
                                                        end if

                                                      end do
                                                     end do
                                                    
        
                                                  end do
                                                 end do
                                               
                                                
                                                ! All the 4 edges have been
                                                ! checked for intersection with
                                                ! the interfacial line
                                                phi_aux(i_p,j_p) = &
                                                SIGN(one,F_tem(i_p,j_p)-0.5)*min( ABS(min_d),ABS(phi_aux(i_p,j_p)) )
                                                       
                                              end if 

                                           end if
                                                     
                                       

                                else
                                end if

                                    end do
                                  end do
                                
                                exit              
                        
                           else
                           end if
                        end do

                         if (flag == 1) then
                         exit
                        else
                        end if

                      end do
                    

                 else
                 end if

              end do
           end do
        


        
          do j = 1,Ny+2,1
             do i = 1,Nx+2,1

                if ((phi_tem(i,j) > 0)) then
                        phi_tt(i,j) = (Kk+1)*d_max
                else if ((phi_tem(i,j) < 0)) then
                        phi_tt(i,j) = -(Kk+1)*d_max
                else
                end if

             end do
          end do
        


        
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1

                flag = 0

                 ! Finding cells where interface exists
                 if ( (F_tem(i,j) > F_low) .and. (F_tem(i,j) < F_up) ) then
                    
                      do j_prime = j-1,j+1,1
                        do i_prime = i-1,i+1,1
                           if ( ( (phi_tem(i,j)*phi_tem(i_prime,j_prime) < 0) &
                              .or. (phi_tem(i,j) == 0) ) .and. ((i-i_prime)*(j-j_prime) <= 1) ) then


                                flag = 1
                                
                                  do j_p = j-Kk,j+Kk,1
                                    do i_p = i-Kk,i+Kk,1

                                    ! Checking that the cells near the interfacial cell (within a band of Kk) indeed are inside the
                                    ! computational domain
                                    if ( (i_p >= 1) .and. (i_p <= Nx+2) .and. (j_p >= 1) .and. &
                                       (j_p <= Ny+2) ) then

                                      phi_tt(i_p,j_p) = phi_aux(i_p,j_p)

                                    else
                                    end if

                                    end do
                                  end do       
                                   
                                exit                      
                        
                           else
                           end if
                        end do

                         if (flag == 1) then
                         exit
                        else
                        end if

                      end do
                    

                 else
                 end if

              end do
           end do
        
        phi_tem = phi_tt

        end subroutine reinitialisePhi

        ! Returns distance of point xd (vector) from the plane in cell (ii,jj)
        function Dis(ii,jj,xd)
                real :: Dis
                integer :: ii,jj
                real, dimension(2) :: xd

                Dis = nor_x(ii,jj)*xd(1) + nor_y(ii,jj)*xd(2) + d_prime(ii,jj)
        end function Dis

        ! Returns Euclidean distance between points p1 and p2 (vectors)
        function Euc(p1,p2)
                real :: Euc
                real, dimension(2) :: p1,p2
                Euc = ( (p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 )**0.5
        end function Euc


        ! To calculate zeta_star (= Fstar .or. phi_star)
        ! zeta_ul = zeta_i-1/2 or zeta_j-1/2 (depending on the sweep direction)
        ! zeta_ur = zeta_i+1/2 or zeta_j+1/2  (zeta = F or phi)
        subroutine calcStar(zeta, zeta_ul, zeta_ur, zeta_star, sweepd)

        implicit none
        integer :: i,j, sweepd
        real, dimension(Nx+2, Ny+2), intent(in) :: zeta, zeta_ul, zeta_ur
        real, dimension(Nx+2, Ny+2), intent(out) :: zeta_star
                
        ! sweepd = 1 for x-sweep, 2 for y-sweep 

        if (sweepd == 1) then

                
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_star(i,j) = ( zeta(i,j) - (del_t/dxe(i))*( zeta_ur(i,j)*u(i,j) - &
                                  zeta_ul(i,j)*u(i-1,j) ) )/( 1 - (del_t/dxe(i))*(u(i,j) - u(i-1,j)) ) 
                                end do
                        end do
              

        else

                
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_star(i,j) = ( zeta(i,j) - (del_t/dye(j))*( zeta_ur(i,j)*v(i,j) - &
                                  zeta_ul(i,j)*v(i,j-1) ) )/( 1 - (del_t/dye(j))*(v(i,j) - v(i,j-1)) )
                                end do
                        end do
              

         end if

        ! Boundary Conditions (Neumann)
        ! x = 0 (Inlet)
        zeta_star(1,:)    = zeta_star(2,:)

        ! x = Lx (Outlet)
        zeta_star(Nx+2,:) = zeta_star(Nx+1,:)

        ! y = 0 (ground)
        zeta_star(:,1)    = zeta_star(:,2)

        ! y = Ly (top)
        zeta_star(:,Ny+2) = zeta_star(:,Ny+1)


        end subroutine calcStar

        ! To calculate zeta_star2 (= Fstar2 .or. phi_star2)
        ! zeta_vb = zeta_i-1/2 or zeta_j-1/2 or zeta_k-1/2 (depending on the sweep direction)
        ! zeta_vt = zeta_i+1/2 or zeta_j+1/2 or zeta_k+1/2 (zeta = F or phi)
        subroutine calc_np1(zeta_star, zeta_vb, zeta_vt, zeta_np1, sweepd)

        implicit none
        integer :: i,j, sweepd
        real, dimension(Nx+2, Ny+2), intent(in) :: zeta_star, zeta_vb, zeta_vt
        real, dimension(Nx+2, Ny+2), intent(out) :: zeta_np1

        ! sweepd = 1 for x-sweep, 2 for y-sweep 

        if (sweepd == 1) then

                
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_np1(i,j) =  zeta_star(i,j)*( 1 + (del_t/dxe(i))*(u(i,j) - u(i-1,j)) ) - &
                                  (del_t/dxe(i))*( zeta_vt(i,j)*u(i,j) - zeta_vb(i,j)*u(i-1,j) )
                                end do
                        end do
                


        else

                
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_np1(i,j) =  zeta_star(i,j)*( 1 + (del_t/dye(j))*(v(i,j) - v(i,j-1)) ) - &
                                  (del_t/dye(j))*( zeta_vt(i,j)*v(i,j) - zeta_vb(i,j)*v(i,j-1) ) 
                                end do
                        end do
               

        end if

        ! Boundary Conditions (Neumann)
        ! x = 0 (Inlet)
        zeta_np1(1,:)    = zeta_np1(2,:)

        ! x = Lx (Outlet)
        zeta_np1(Nx+2,:) = zeta_np1(Nx+1,:)

        ! y = 0 (ground)
        zeta_np1(:,1)    = zeta_np1(:,2)

        ! y = Ly (top)
        zeta_np1(:,Ny+2) = zeta_np1(:,Ny+1)


        end subroutine calc_np1



        ! To calculate phi on i-1/2 and i+1/2 faces
        ! phi_sweep = phi (x-y sweep)/phi_star (y-x sweep)
        subroutine phiFaceCenter_xSweep(phi_sweep)

        implicit none
        integer :: i,j
        real, dimension(Nx+2, Ny+2) :: phi_sweep

        
                do j = 2,Ny+1,1
                        do i = 1,Nx+1,1

                                ! right face
                                if (u(i,j) >= 0) then

                                        if (i == 1) then
                                                phi_ur(i,j) = phi_sweep(i,j) + (dxe(i)/2)*(1 - &
                                                u(i,j)*(del_t/dxe(i)))*((phi_sweep(i+1,j) - &
                                                phi_sweep(i,j))/(xc(i+1) - xc(i)))
                                        else
                                                phi_ur(i,j) = phi_sweep(i,j) + (dxe(i)/2)*(1 - &
                                                u(i,j)*(del_t/dxe(i)))*((phi_sweep(i+1,j) - &
                                                phi_sweep(i-1,j))/(xc(i+1) - xc(i-1)))
                                        end if

                                else 

                                        if (i == Nx+1) then
                                                phi_ur(i,j) = phi_sweep(i,j) + (dxe(i)/2)*(1 + &
                                                u(i,j)*(del_t/dxe(i)))*((phi_sweep(i+1,j) - &
                                                phi_sweep(i,j))/(xc(i+1) - xc(i)))
                                        else
                                                phi_ur(i,j) = phi_sweep(i,j) + (dxe(i)/2)*(1 + &
                                                u(i,j)*(del_t/dxe(i)))*((phi_sweep(i+2,j) - &
                                                phi_sweep(i,j))/(xc(i+2) - xc(i)))
                                        end if

                                end if
                        
                                ! left face
                                phi_ul(i+1,j) = phi_ur(i,j)                        
                        end do
                end do
       

        end subroutine phiFaceCenter_xSweep

        ! To calculate phi on j-1/2 and j+1/2 faces
        ! phi_sweep = phi (y-x sweep)/phi_star (x-y zweep)
        subroutine phiFaceCenter_ySweep(phi_sweep)

        implicit none
        integer :: i,j
        real, dimension(Nx+2, Ny+2) :: phi_sweep

        
                do j = 1,Ny+1,1
                        do i = 2,Nx+1,1

                                ! top face
                                if (v(i,j) >= 0) then

                                        if (j == 1) then
                                                phi_vt(i,j) = phi_sweep(i,j) + (dye(j)/2)*(1 - &
                                                v(i,j)*(del_t/dye(j)))*((phi_sweep(i,j+1) - &
                                                phi_sweep(i,j))/(yc(j+1) - yc(j)))
                                        else
                                                phi_vt(i,j) = phi_sweep(i,j) + (dye(j)/2)*(1 - &
                                                v(i,j)*(del_t/dye(j)))*((phi_sweep(i,j+1) - &
                                                phi_sweep(i,j-1))/(yc(j+1) - yc(j-1)))
                                        end if

                                else 

                                        if (j == Ny+1) then
                                                phi_vt(i,j) = phi_sweep(i,j) + (dye(j)/2)*(1 + &
                                                v(i,j)*(del_t/dye(j)))*((phi_sweep(i,j+1) - &
                                                phi_sweep(i,j))/(yc(j+1) - yc(j)))
                                        else
                                                phi_vt(i,j) = phi_sweep(i,j) + (dye(j)/2)*(1 + &
                                                v(i,j)*(del_t/dye(j)))*((phi_sweep(i,j+2) - &
                                                phi_sweep(i,j))/(yc(j+2) - yc(j)))
                                        end if

                                end if

                                ! bottom face
                                phi_vb(i,j+1) = phi_vt(i,j)

                        end do
                end do
      

        end subroutine phiFaceCenter_ySweep


        ! To solve Ax = b where A = n*n matrix
        ! Note that this subroutine modifies the input matrices A and b as well
        subroutine GaussElimination(A,x,b,n)

        implicit none
        integer :: n
        real, dimension(n)               :: x
        real, dimension(n)               :: b
        real, dimension(n,n)             :: A
        real, dimension(n)               :: ss
        integer :: ii,jj,k_k,p
        real :: big, dummy, factor, summ

        do ii = 1,n,1
                ss(ii) = ABS(A(ii,1))
                do jj = 2,n,1
                        if (ABS(A(ii,jj)) > ss(ii)) then
                                ss(ii) = ABS(A(ii,jj))
                        else
                        end if
                end do
        end do

        ! Forward Elimination
        do k_k = 1,n-1,1

                ! Partial Pivoting
                p = k_k
                big = ABS(A(k_k,k_k)/ss(k_k))
                do ii = k_k+1,n,1
                        dummy = ABS(A(ii,k_k)/ss(ii))
                        if (dummy > big) then
                                big = dummy
                                p = ii
                        else
                        end if
                end do
                
                if (p /= k_k) then
                        do jj = k_k,n,1
                                dummy = A(p,jj)
                                A(p,jj) = A(k_k,jj)
                                A(k_k,jj) = dummy
                        end do
                        dummy = b(p)
                        b(p) = b(k_k)
                        b(k_k) = dummy
                        dummy = ss(p)
                        ss(p) = ss(k_k)
                        ss(k_k) = dummy
                else
                end if
                
                do ii = k_k+1,n,1
                        factor = A(ii,k_k)/A(k_k,k_k)
                        do jj = k_k+1,n,1
                                A(ii,jj) = A(ii,jj) - factor*A(k_k,jj)
                        end do
                        b(ii) = b(ii) - factor*b(k_k)
                end do

        end do

        ! Back Substitution
        x(n) = b(n)/A(n,n)
        do ii = n-1,1,-1
                summ = 0
                do jj = ii+1,n,1
                        summ = summ + A(ii,jj)*x(jj) 
                end do
                x(ii) = (b(ii) - summ)/A(ii,ii)
        end do        

        end subroutine GaussElimination

        

        subroutine delta_t(delt)
        
                implicit none
                integer :: ikea, jkea
                real :: delt, delt_x, delt_y, delx, dely, umax, vmax
                umax = 0
                vmax = 0

                delx = minval(dxe(2:Nx+2))
                dely = minval(dye(2:Ny+2))

                
                  do jkea = 2,Ny+1,1
                    do ikea = 2,Nx+1,1
                        
                        if (umax < ABS(u(ikea,jkea))) then
                                umax = ABS(u(ikea,jkea))
                        else
                        end if

                        if (vmax < ABS(v(ikea,jkea))) then
                                vmax = ABS(v(ikea,jkea))
                        else
                        end if

                    end do
                  end do
                

                delt_x = (0.25*delx)/umax
                delt_y = (0.25*dely)/vmax

                if ((umax < (0.25*delx)) .and. (vmax < 0.25*dely)) then
                        delt = 0.01
                else        
                        delt = min(delt_x, delt_y)
                end if

        end subroutine delta_t



end program main
