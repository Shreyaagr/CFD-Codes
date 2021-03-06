! ME634A end-sem assignment : Three-dimensional stretching of a sphere in a shear velocity field using CLSVOF method.
! Code by Shreya Agrawal (160662) in Fortran 95

program main
        implicit none
        
        real :: Lx, Ly, Lz, del_t, tol, time, T, radius, center_x, center_y, center_z, dVf, resNR, tolNR, eps, min_d
	real :: flag, F_low, F_up
        real, dimension(:,:,:), allocatable :: u, v, w, F, Fstar, Fstar2, F_np1, phi, phi_star, &
                                                phi_star2, phi_np1, F_ul, F_ur, F_vb, F_vt, F_wn, F_ws, phi_aux, &
                                                phi_ul, phi_ur, phi_vb, phi_vt, phi_wn, phi_ws, nor_x, nor_y, nor_z, d, s, d_prime 
        integer :: Nx, Ny, Nz, i, j, k, t_its, temp, initialF, itsNR, k_pp, j_pp, i_pp
        real, dimension(:), allocatable :: xe, xc, ye, yc, ze, zc, dxe, dxc, dye, dyc, dze, dzc  ! grid edge and center coordinates for each grid-cell
        real, parameter :: pi = 3.1415927
        real, parameter :: one = 1.0


        ! Time-period for velocity
        T = 3

        ! Compuatational domain
        Lx = 1
        Ly = 1
        Lz = 1

        ! Cylinder radius and center
        radius   = 0.15
        center_x = 0.35
        center_y = 0.35
        center_z = 0.35
        
        ! Number of internal nodes in x,y,z directions
        Nx = 64
        Ny = 64
        Nz = 64

        ! Current time and iteration no. in time
        time  = 0
        t_its = 1

        min_d = 10000
	flag = 0

	F_low = 1.0e-6
	F_up = 1.0 - 1.0e-6

        temp = 0
        initialF = 0 ! Make =1 when initialising F (t = 0), and 0 at all other times 

        resNR = 1       ! Residual for NR iterations
        tolNR = 1.0e-5  ! tolerance for NR iterations
        itsNR = 0       ! NR iteration counter
        
        ! Thickness of the interface = epsilon
        eps = 1.75*(Lx/Nx)

        ! Matrix size allocation
        allocate(u(Nx+2, Ny+2, Nz+2))   ! x-velocity
        allocate(v(Nx+2, Ny+2, Nz+2))   ! y-velocity
        allocate(w(Nx+2, Ny+2, Nz+2))   ! z-velocity
        allocate(F(Nx+2, Ny+2, Nz+2))   ! VOF function
        allocate(phi(Nx+2, Ny+2, Nz+2))   ! LS function

        allocate(Fstar(Nx+2, Ny+2, Nz+2))   
        allocate(phi_star(Nx+2, Ny+2, Nz+2))

        allocate(Fstar2(Nx+2, Ny+2, Nz+2))
        allocate(phi_star2(Nx+2, Ny+2, Nz+2))

        allocate(F_np1(Nx+2, Ny+2, Nz+2))   ! VOF function at next time-step
        allocate(phi_np1(Nx+2, Ny+2, Nz+2)) ! LS function at next time-step

        allocate(phi_aux(Nx+2, Ny+2, Nz+2))

        allocate(F_ul(Nx+2, Ny+2, Nz+2))    ! VOF function at the cell edges
        allocate(F_ur(Nx+2, Ny+2, Nz+2))
        allocate(F_vb(Nx+2, Ny+2, Nz+2))
        allocate(F_vt(Nx+2, Ny+2, Nz+2))
        allocate(F_wn(Nx+2, Ny+2, Nz+2))
        allocate(F_ws(Nx+2, Ny+2, Nz+2))

        allocate(phi_ul(Nx+2, Ny+2, Nz+2))    ! LS function at the cell edges
        allocate(phi_ur(Nx+2, Ny+2, Nz+2))
        allocate(phi_vb(Nx+2, Ny+2, Nz+2))
        allocate(phi_vt(Nx+2, Ny+2, Nz+2))
        allocate(phi_wn(Nx+2, Ny+2, Nz+2))
        allocate(phi_ws(Nx+2, Ny+2, Nz+2))

        allocate(nor_x(Nx+2, Ny+2, Nz+2))
        allocate(nor_y(Nx+2, Ny+2, Nz+2))
        allocate(nor_z(Nx+2, Ny+2, Nz+2))
        allocate(d(Nx+2, Ny+2, Nz+2))
        allocate(s(Nx+2, Ny+2, Nz+2))
        allocate(d_prime(Nx+2, Ny+2, Nz+2))

        allocate(xe(Nx+2))
        allocate(xc(Nx+2))
        allocate(ye(Ny+2))
        allocate(yc(Ny+2))
        allocate(ze(Nz+2))
        allocate(zc(Nz+2))

        allocate(dxe(Nx+2))
        allocate(dxc(Nx+2))
        allocate(dye(Ny+2))
        allocate(dyc(Ny+2))
        allocate(dze(Nz+2))
        allocate(dzc(Nz+2))

        ! Initialisation
        u = 0
        v = 0
        w = 0
        F = 0
        phi = 0

        Fstar  = 0
        Fstar2 = 0
        F_np1  = 0 

        phi_star  = 0
        phi_star2 = 0
        phi_np1   = 0

        phi_aux = 10000

        F_ul = 0
        F_ur = 0
        F_vb = 0
        F_vt = 0
        F_wn = 0
        F_ws = 0

        phi_ul = 0
        phi_ur = 0
        phi_vb = 0
        phi_vt = 0
        phi_wn = 0
        phi_ws = 0

        nor_x = 0
        nor_y = 0
        nor_z = 0
        d     = 0
        s     = 0
        d_prime = 0

        dxe = 0
        dxc = 0
        dye = 0
        dyc = 0
        dze = 0
        dzc = 0

        ! Presently the grid is uniform, but can be made non-uniform through the following assignments
        do i = 1,Nx+2,1
                xe(i) = (i-1)*(Lx/Nx)
                xc(i) = (i-2+0.5)*(Lx/Nx)
        end do
        
        do j = 1,Ny+2,1
                ye(j) = (j-1)*(Ly/Ny)
                yc(j) = (j-2+0.5)*(Ly/Ny)
        end do

        do k = 1,Nz+2,1
                ze(k) = (k-1)*(Lz/Nz)
                zc(k) = (k-2+0.5)*(Lz/Nz)
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

        do k = 2,Nz+2,1
               dze(k) = ze(k) - ze(k-1)
        end do

        do k = 1,Nz+1,1
               dzc(k) = zc(k+1) - zc(k)
        end do

        ! Ghost Cells grid-size
        dxe(1)    = dxe(2)
        dxc(Nx+2) = dxc(Nx+1)
        dye(1)    = dye(2)
        dyc(Ny+2) = dyc(Ny+1)
        dze(1)    = dze(2)
        dzc(Nz+2) = dzc(Nz+1)

                ! velocity calculation for time = 0 
                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                        u(i,j,k) =  2*((SIN(pi*xe(i)))**2)*SIN(2*pi*yc(j))*SIN(2*pi*zc(k))*COS((pi*time)/T) 
                                        v(i,j,k) =  -((SIN(pi*ye(j)))**2)*SIN(2*pi*xc(i))*SIN(2*pi*zc(k))*COS((pi*time)/T)
                                        w(i,j,k) =  -((SIN(pi*ze(k)))**2)*SIN(2*pi*xc(i))*SIN(2*pi*yc(j))*COS((pi*time)/T)
                                end do
                        end do
                end do

        ! velocity BCs
        ! x = 0 (Inlet)
        u(1,:,:) = 0
        v(1,:,:) = -v(2,:,:)
        w(1,:,:) = -w(2,:,:)

        ! x = Lx (Outlet)
        u(Nx+1,:,:) = 0
        v(Nx+2,:,:) = -v(Nx+1,:,:)
        w(Nx+2,:,:) = -w(Nx+1,:,:)

        ! y = 0 (ground)
        v(:,1,:) = 0
        u(:,1,:) = -u(:,2,:)
        w(:,1,:) = -w(:,2,:)

        ! y = Ly (top)
        v(:,Ny+1,:) = 0
        u(:,Ny+2,:) = -u(:,Ny+1,:)
        w(:,Ny+2,:) = -w(:,Ny+1,:)

        ! z = 0 (front)
        w(:,:,1) = 0
        u(:,:,1) = -u(:,:,2)
        v(:,:,1) = -v(:,:,2)

        ! z = Lz (back)
        w(:,:,Nz+1) = 0
        u(:,:,Nz+2) = -u(:,:,Nz+1)
        v(:,:,Nz+2) = -v(:,:,Nz+1)


        
        ! LS function at t = 0 
        ! Pg.978 of Ref7
        ! Check this - phi<0 inside circle and phi>0 outside circle
        do k = 2,Nz+1,1
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1
                phi(i,j,k) = -radius + SQRT((xc(i) - center_x)**2 + (yc(j) - center_y)**2 + &
                                                (zc(k) - center_z)**2)
              end do
           end do
        end do
        

        ! Boundary Conditions for LS function (Neumann)
        ! ******* CONFIRM ***********

        ! x = 0 (Inlet)
        phi(1,:,:) = phi(2,:,:)

        ! x = Lx (Outlet)
        phi(Nx+2,:,:) = phi(Nx+1,:,:)

        ! y = 0 (ground)
        phi(:,1,:) = phi(:,2,:)

        ! y = Ly (top)
        phi(:,Ny+2,:) = phi(:,Ny+1,:)

        ! z = 0 (front)
        phi(:,:,1) = phi(:,:,2)

        ! z = Lz (back)
        phi(:,:,Nz+2) = phi(:,:,Nz+1)


        ! Calculating VOF function at t = 0 
        initialF = 1
        call reconstruct(F, phi)

        ! Initialising F in all cells
        do k = 2,Nz+1,1
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1

                if (phi(i,j,k) <= 0) then
                   F(i,j,k) = 0
                else
                   F(i,j,k) = 1
                end if

              end do
           end do
        end do

        ! Correcting F in interfacial cells
        do k = 2,Nz+1,1
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1

		flag = 0

                ! Finding cells where interface exists
                do k_pp = k-1,k+1,1
                  do j_pp = j-1,j+1,1
                    do i_pp = i-1,i+1,1
                       if ( ( (phi(i,j,k)*phi(i_pp,j_pp,k_pp) < 0) &
                          .or. (phi(i,j,k) == 0) ) .and. (ABS(phi(i,j,k)) <= eps) &
                        .and. ((i-i_pp)*(j-j_pp)*(k-k_pp) <= 1) ) then

                   dVf = del_Vf(nor_x(i,j,k),nor_y(i,j,k),nor_z(i,j,k),s(i,j,k),dxe(i),dye(j),dze(k), phi(i,j,k))
                   F(i,j,k) = dVf/(dxe(i)*dye(j)*dze(k)) 

			flag = 1

                       else
                       end if

                    end do
			
			 if (flag == 1) then
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
        end do

       

        ! Boundary Conditions for VOF function (Neumann)
        ! ******* CONFIRM ***********

        ! x = 0 (Inlet)
        F(1,:,:)   = F(2,:,:)

        ! x = Lx (Outlet)
        F(Nx+2,:,:)   = F(Nx+1,:,:)

        ! y = 0 (ground)
        F(:,1,:)   = F(:,2,:)

        ! y = Ly (top)
        F(:,Ny+2,:)   = F(:,Ny+1,:)

        ! z = 0 (front)
        F(:,:,1)   = F(:,:,2)

        ! z = Lz (back) 
        F(:,:,Nz+2)   = F(:,:,Nz+1)

                ! Write to output files
                ! LS function
                open(unit = 100, file = 'phi_N64_t0.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(100,*) phi(i,j,k)
                                end do
                        end do
                end do
                close(100)

                ! VOF function
                open(unit = 200, file = 'F_N64_t0.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(200,*) F(i,j,k)
                                end do
                        end do
                end do
                close(200)

        

        initialF = 0 ! VOF initialisation is complete


        open(unit = 1, file = 'fort.200')

        ! Solving in time
        do while ( time <= T )

                ! del_t calculation
                call delta_t(del_t)

                ! update present time
                time = time + del_t

                write(*,*) 'time =',time

                ! velocity calculation
                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                        u(i,j,k) = 2*((SIN(pi*xe(i)))**2)*SIN(2*pi*yc(j))*SIN(2*pi*zc(k))*COS((pi*time)/T)
                                        v(i,j,k) = -((SIN(pi*ye(j)))**2)*SIN(2*pi*xc(i))*SIN(2*pi*zc(k))*COS((pi*time)/T)
                                        w(i,j,k) = -((SIN(pi*ze(k)))**2)*SIN(2*pi*xc(i))*SIN(2*pi*yc(j))*COS((pi*time)/T)         
                                end do
                        end do
                end do                

        ! velocity BCs
        ! x = 0 (Inlet)
        u(1,:,:) = 0
        v(1,:,:) = -v(2,:,:)
        w(1,:,:) = -w(2,:,:)

        ! x = Lx (Outlet)
        u(Nx+1,:,:) = 0
        v(Nx+2,:,:) = -v(Nx+1,:,:)
        w(Nx+2,:,:) = -w(Nx+1,:,:)

        ! y = 0 (ground)
        v(:,1,:) = 0
        u(:,1,:) = -u(:,2,:)
        w(:,1,:) = -w(:,2,:)

        ! y = Ly (top)
        v(:,Ny+1,:) = 0
        u(:,Ny+2,:) = -u(:,Ny+1,:)
        w(:,Ny+2,:) = -w(:,Ny+1,:)

        ! z = 0 (front)
        w(:,:,1) = 0
        u(:,:,1) = -u(:,:,2)
        v(:,:,1) = -v(:,:,2)

        ! z = Lz (back)
        w(:,:,Nz+1) = 0
        u(:,:,Nz+2) = -u(:,:,Nz+1)
        v(:,:,Nz+2) = -v(:,:,Nz+1) 
         
   
                  
        if ( MOD(t_its - 3*temp,3) == 0 ) then
        ! z-x-y sweep
                write(*,*) 'Going z-x-y!'

                ! Calculating phi_wn, phi_ws   **** CONFIRM THIS ****
                call phiFaceCenter_zSweep(phi)
                call reconstruct(F,phi)
                call FfaceCenter_zSweep(F,phi)
                call calcStar(phi,phi_wn,phi_ws,phi_star,3)
                call calcStar(F,F_wn,F_ws,Fstar,3)

                call phiFaceCenter_xSweep(phi_star)
                call reconstruct(Fstar,phi_star)
                call FfaceCenter_xSweep(Fstar,phi_star)
                call calcStar2(phi_star,phi_ul,phi_ur,phi_star2,1)
                call calcStar2(Fstar,F_ul,F_ur,Fstar2,1)

                call phiFaceCenter_ySweep(phi_star2)
                call reconstruct(Fstar2,phi_star2)
                call FfaceCenter_ySweep(Fstar2, phi_star2)
                call calc_np1(phi_star,phi_star2,phi_vb,phi_vt,phi_np1,2)
                call calc_np1(Fstar,Fstar2,F_vb,F_vt,F_np1,2)


                ! Truncating VOF Function
                do k = 1,Nz+2,1
                   do j = 1,Ny+2,1
                      do i = 1,Nx+2,1

                        if ( (phi_np1(i,j,k) < -eps) .or. (F_np1(i,j,k) < F_low) ) then
                                F_np1(i,j,k) = 0
                        else if ( (phi_np1(i,j,k) > eps) .or. (F_np1(i,j,k) > F_up) ) then
                                F_np1(i,j,k) = 1
                        else
                        end if

                      end do
                   end do
                end do

                call reconstruct(F_np1,phi_np1)
                phi_aux = 10000
                call reinitialisePhi(F_np1,phi_np1)


        else if ( MOD(t_its - 3*temp,2) == 0 ) then
        ! y-z-x sweep
                write(*,*) 'Going y-z-x!'

                call phiFaceCenter_ySweep(phi)
                call reconstruct(F,phi)
                call FfaceCenter_ySweep(F, phi)
                call calcStar(phi,phi_vb,phi_vt,phi_star,2)
                call calcStar(F,F_vb,F_vt,Fstar,2)

                call phiFaceCenter_zSweep(phi_star)
                call reconstruct(Fstar,phi_star)
                call FfaceCenter_zSweep(Fstar,phi_star)
                call calcStar2(phi_star,phi_wn,phi_ws,phi_star2,3)
                call calcStar2(Fstar,F_wn,F_ws,Fstar2,3)

                call phiFaceCenter_xSweep(phi_star2)
                call reconstruct(Fstar2,phi_star2)
                call FfaceCenter_xSweep(Fstar2,phi_star2)
                call calc_np1(phi_star,phi_star2,phi_ul,phi_ur,phi_np1,1)
                call calc_np1(Fstar,Fstar2,F_ul,F_ur,F_np1,1)

                ! Truncating VOF Function
                do k = 1,Nz+2,1
                   do j = 1,Ny+2,1
                      do i = 1,Nx+2,1

                        if ( (phi_np1(i,j,k) < -eps) .or. (F_np1(i,j,k) < F_low) ) then
                                F_np1(i,j,k) = 0
                        else if ( (phi_np1(i,j,k) > eps) .or. (F_np1(i,j,k) > F_up) ) then
                                F_np1(i,j,k) = 1
                        else
                        end if

                      end do
                   end do
                end do

                call reconstruct(F_np1,phi_np1)
                phi_aux = 10000
                call reinitialisePhi(F_np1,phi_np1)


        else
        ! x-y-z sweep
                write(*,*) 'Going x-y-z!'
                
                call phiFaceCenter_xSweep(phi)
                call reconstruct(F,phi)
                call FfaceCenter_xSweep(F,phi)
                call calcStar(phi,phi_ul,phi_ur,phi_star,1)
                call calcStar(F,F_ul,F_ur,Fstar,1)

                call phiFaceCenter_ySweep(phi_star)
                call reconstruct(Fstar,phi_star)
                call FfaceCenter_ySweep(Fstar,phi_star)
                call calcStar2(phi_star,phi_vb,phi_vt,phi_star2,2)
                call calcStar2(Fstar,F_vb,F_vt,Fstar2,2)

                call phiFaceCenter_zSweep(phi_star2)
                call reconstruct(Fstar2,phi_star2)
                call FfaceCenter_zSweep(Fstar2,phi_star2)
                call calc_np1(phi_star,phi_star2,phi_wn,phi_ws,phi_np1,3)
                call calc_np1(Fstar,Fstar2,F_wn,F_ws,F_np1,3)

                ! Truncating VOF Function
                do k = 1,Nz+2,1
                   do j = 1,Ny+2,1
                      do i = 1,Nx+2,1
                        
                        if ( (phi_np1(i,j,k) < -eps) .or. (F_np1(i,j,k) < F_low) ) then
                                F_np1(i,j,k) = 0
                        else if ( (phi_np1(i,j,k) > eps) .or. (F_np1(i,j,k) > F_up) ) then
                                F_np1(i,j,k) = 1
                        else
                        end if

                      end do
                   end do
                end do

                call reconstruct(F_np1,phi_np1)
                phi_aux = 10000
                call reinitialisePhi(F_np1,phi_np1)

        end if
          
        if ( MOD(t_its,3) == 0 ) then
                temp = temp + 1
        else
        end if
  
                phi = phi_np1
                F   = F_np1

                t_its = t_its + 1

                ! Write to output files

           if ( (time >= 0.47) .and. (time <= 0.49) ) then

                ! LS function
                open(unit = 300, file = 'phi_N64_t0pt48.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(300,*) phi(i,j,k)
                                end do
                        end do
                end do
                close(300)

                ! VOF function
                open(unit = 400, file = 'F_N64_t0pt48.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(400,*) F(i,j,k)
                                end do
                        end do
                end do
                close(400)

           else if ( (time >= 0.97) .and. (time <= 0.99) ) then

                ! Write to output files
                ! LS function
                open(unit = 500, file = 'phi_N64_t0pt98.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(500,*) phi(i,j,k)
                                end do
                        end do
                end do
                close(500)

                ! VOF function
                open(unit = 600, file = 'F_N64_t0pt98.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(600,*) F(i,j,k)
                                end do
                        end do
                end do
                close(600)

           else if ( (time >= 1.47) .and. (time <= 1.49) ) then

                ! Write to output files
                ! LS function
                open(unit = 700, file = 'phi_N64_t1pt49.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(700,*) phi(i,j,k)
                                end do
                        end do
                end do
                close(700)

                ! VOF function
                open(unit = 800, file = 'F_N64_t1pt49.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(800,*) F(i,j,k)
                                end do
                        end do
                end do
                close(800)

           else if ( (time >= 1.97) .and. (time <= 1.99) ) then

                ! Write to output files
                ! LS function
                open(unit = 900, file = 'phi_N64_t1pt98.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(900,*) phi(i,j,k)
                                end do
                        end do
                end do
                close(900)

                ! VOF function
                open(unit = 1000, file = 'F_N64_t1pt98.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(1000,*) F(i,j,k)
                                end do
                        end do
                end do
                close(1000)     

           else if ( (time >= 2.47) .and. (time <= 2.49) ) then

                ! Write to output files
                ! LS function
                open(unit = 1100, file = 'phi_N64_t2pt48.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(1100,*) phi(i,j,k)
                                end do
                        end do
                end do
                close(1100)

                ! VOF function
                open(unit = 1200, file = 'F_N64_t2pt48.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(1200,*) F(i,j,k)
                                end do
                        end do
                end do
                close(1200)

           else
           end if


        end do


                ! Write to output files
                ! LS function
                open(unit = 1300, file = 'phi_N64_t3.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(1300,*) phi(i,j,k)
                                end do
                        end do
                end do
                close(1300)

                ! VOF function
                open(unit = 1400, file = 'F_N64_t3.txt')
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                write(1400,*) F(i,j,k)
                                end do
                        end do
                end do
                close(1400)


                close(1)

contains

        ! To reconstruct the interface using phi and F
        ! Gives nor_x, nor_y, nor_z, d (without mass conservation) and s (ensuring mass conservation)
        ! F_temp = F/Fstar/Fstar2/F_np1, phi_temp = phi/phi_star/phi_star2/phi_np1
        subroutine reconstruct(F_temp, phi_temp)

        implicit none
        integer :: i, j, k, i_prime, j_prime, k_prime, i_p, j_p, k_p
        real :: tt, Fc, dx1, dy1, dz1, t1, t2, t3, sc, s1, sm, w8, Hp, Xt, Yt, Zt 
        real, dimension(Nx+2, Ny+2, Nz+2), intent(in) :: F_temp, phi_temp
        real, dimension(4,4) :: A
        real, dimension(4)   :: x, b        
        real, dimension(3,3) :: A3
        real, dimension(3)   :: x3, b3
        

        tt = 0 ! use tt to normalise the plane 
        
        do k = 2,Nz+1,1
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1
                            
                 ! Assign zero normals to all cells and update for only
                 ! interfacial cells
                 nor_x(i,j,k) = 0
                 nor_y(i,j,k) = 0
                 nor_z(i,j,k) = 0
                 d(i,j,k)     = 0
                 s(i,j,k)     = 0

		flag = 0

                 ! Finding cells where interface exists        
                 if ( ((F_temp(i,j,k) > F_low) .and. (F_temp(i,j,k) < F_up)) .or. (initialF == 1)  ) then
                    do k_prime = k-1,k+1,1
                      do j_prime = j-1,j+1,1
                        do i_prime = i-1,i+1,1
                           if ( ( (phi_temp(i,j,k)*phi_temp(i_prime,j_prime,k_prime) < 0) .or. &
                                (phi_temp(i,j,k) == 0) ) .and. ( ABS(phi_temp(i,j,k)) <= eps )  &
                                .and. ((i-i_prime)*(j-j_prime)*(k-k_prime) <= 1) ) then

                		  flag = 1                                                                                                
                 
		                  ! Reconstruction using phi 

                                   A = 0
                                   b = 0
                                   x = 0
                                   A3 = 0
                                   b3 = 0
                                   x3 = 0
        
                                   do k_p = k-1,k+1,1
                                      do j_p = j-1,j+1,1
                                         do i_p = i-1,i+1,1

                                            Xt = xc(i_p) - xc(i)
                                            Yt = yc(j_p) - yc(j)
                                            Zt = zc(k_p) - zc(k)

                                            if ((i_p == i) .and. (j_p == j) .and. (k_p == k)) then
                                                w8 = 52
                                            else
                                                w8 = 1
                                            end if
                                                
                                            if ( ABS(phi_temp(i_p,j_p,k_p)) <= eps ) then
                                                Hp = (1/(2*eps))*(1+COS((pi*phi_temp(i_p,j_p,k_p))/eps))
                                            else
                                                Hp = 0
                                            end if    
                                        

                                            A(1,1) = A(1,1) + w8*Hp*Xt*Xt
                                            A(2,1) = A(2,1) + w8*Hp*Xt*Yt
                                            A(3,1) = A(3,1) + w8*Hp*Xt*Zt
                                            A(4,1) = A(4,1) + w8*Hp*Xt

                                            A(1,2) = A(1,2) + w8*Hp*Xt*Yt
                                            A(2,2) = A(2,2) + w8*Hp*Yt*Yt
                                            A(3,2) = A(3,2) + w8*Hp*Yt*Zt
                                            A(4,2) = A(4,2) + w8*Hp*Yt
                                            
                                            A(1,3) = A(1,3) + w8*Hp*Xt*Zt
                                            A(2,3) = A(2,3) + w8*Hp*Yt*Zt
                                            A(3,3) = A(3,3) + w8*Hp*Zt*Zt
                                            A(4,3) = A(4,3) + w8*Hp*Zt

                                            A(1,4) = A(1,4) + w8*Hp*Xt
                                            A(2,4) = A(2,4) + w8*Hp*Yt
                                            A(3,4) = A(3,4) + w8*Hp*Zt
                                            A(4,4) = A(4,4) + w8*Hp

                                            b(1) = b(1) + w8*Hp*phi_temp(i_p,j_p,k_p)*Xt
                                            b(2) = b(2) + w8*Hp*phi_temp(i_p,j_p,k_p)*Yt
                                            b(3) = b(3) + w8*Hp*phi_temp(i_p,j_p,k_p)*Zt
                                            b(4) = b(4) + w8*Hp*phi_temp(i_p,j_p,k_p)
                                                
                                         end do
                                      end do
                                   end do       

                                                                           
                                   ! Use Gauss Elimination to find the plane normals
                                        call GaussElimination(A,x,b,4)

                                   ! Normalise x
                                   tt = ( x(1)**2 + x(2)**2 + x(3)**2 )**0.5
                                   nor_x(i,j,k) = x(1)/tt
                                   nor_y(i,j,k) = x(2)/tt
                                   nor_z(i,j,k) = x(3)/tt
                                   d(i,j,k)     = x(4)/tt
                                   
                                   ! Calculate s using d - Pg.11 of CLSVOF doc       
                                   s(i,j,k) = d(i,j,k) + 0.5*( ABS(nor_x(i,j,k))*dxe(i) + &
                                              ABS(nor_y(i,j,k))*dye(j) + ABS(nor_z(i,j,k))*dze(k) )     

                                   ! Correcting s using F
                                   if (initialF == 0) then
                                         
                                       t1 = dxe(i)*ABS(nor_x(i,j,k))
                                       t2 = dye(j)*ABS(nor_y(i,j,k))
                                       t3 = dze(k)*ABS(nor_z(i,j,k))

                                       ! Ensuring that dx1t >= dy1t >= dz1t
                                       dx1 = max(t1,t2,t3)
                                       dz1 = min(t1,t2,t3)
                                       dy1 = t1 + t2 + t3 - dx1 - dz1

                                       sm = dx1 + dy1 + dz1
                                       sc = min(s(i,j,k),sm-s(i,j,k)) 
                                       Fc = min(F_temp(i,j,k),1-F_temp(i,j,k)) 
                                        
                                       sc = (6*Fc*dx1*dy1*dz1)**(1.0/3.0)
                                       s1 = sc
 
                                       if (sc < dz1) then
                                         if (F_temp(i,j,k) <= 0.5) then
                                             s(i,j,k) = sc
                                         else 
                                             s(i,j,k) = sm - sc
                                         end if
                                       else 
                                         if (2*Fc*dx1*dy1 - ((dz1**2)/12) >= 0) then
                                            sc = 0.5*dz1 + (2*Fc*dx1*dy1 - ((dz1**2)/12))**0.5

                                            s1 = sc
                                         else
                                         end if
                                        
                                         if ((sc < dy1) .and. (2*Fc*dx1*dy1 - ((dz1**2)/12) >= 0)) then
                                           if (F_temp(i,j,k) <= 0.5) then
                                             s(i,j,k) = sc
                                           else
                                             s(i,j,k) = sm - sc
                                           end if

                                         else
                                           if (dx1 >= dy1+dz1) then
                                             sc = Fc*dx1 + (dy1+dz1)/2
                                             if (sc >= dy1+dz1) then
                                                if (F_temp(i,j,k) <= 0.5) then
                                                   s(i,j,k) = sc
                                                else
                                                   s(i,j,k) = sm - sc
                                                end if
                                             else
                                                ! NR iterations
                                                sc = s1 ! Initial guess

                                                resNR = 1
                                                itsNR = 0
                                                call NewtonR(sc,Fc,dx1,dy1,dz1)                                         

                                                if (F_temp(i,j,k) <= 0.5) then
                                                   s(i,j,k) = sc
                                                else
                                                   s(i,j,k) = sm - sc
                                                end if

                                             end if
                                           else
                                                ! NR iterations
                                                sc = s1 ! initial guess

                                                resNR = 1
                                                itsNR = 0
                                                call NewtonR(sc,Fc,dx1,dy1,dz1)
                                                
                                                if (F_temp(i,j,k) <= 0.5) then
                                                   s(i,j,k) = sc
                                                else
                                                   s(i,j,k) = sm - sc
                                                end if

                                           end if
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

			 if (flag == 1) then
                         exit
                        else
                        end if

                    end do
                                     
                 else
                 end if

              end do
           end do
        end do        

        end subroutine reconstruct

        ! Newton Raphson iterations
        subroutine NewtonR(sc_t,Fc_t,dx1_t,dy1_t,dz1_t)
        
        implicit none       
        real, intent(in) :: Fc_t, dx1_t, dy1_t, dz1_t 
        real, intent(out):: sc_t
        real :: sc_tnew

        resNR = 1
        itsNR = 0  ! No. of iterations
        
        do while ((resNR > tolNR) .and. (itsNR < 50))

        if (sc_t < dx1_t) then

                sc_tnew = sc_t - (6*Fc_t*dx1_t*dy1_t*dz1_t - sc_t**3 + (sc_t - dz1_t)**3 + &
                (sc_t - dy1_t)**3)/(-3*(sc_t**2) + 3*((sc_t - dz1_t)**2) + 3*((sc_t - dy1_t)**2))                

                resNR = ABS( 6*Fc_t*dx1_t*dy1_t*dz1_t - sc_tnew**3 + (sc_tnew - dz1_t)**3 + &
                        (sc_tnew - dy1_t)**3 )

        else

                sc_tnew = sc_t - (6*Fc_t*dx1_t*dy1_t*dz1_t - sc_t**3 + (sc_t - dz1_t)**3 + &
                (sc_t - dy1_t)**3 + (sc_t - dx1_t)**3)/(-3*(sc_t**2) + 3*((sc_t - dz1_t)**2) + &
                3*((sc_t - dy1_t)**2) + 3*((sc_t - dx1_t)**2))

        
                resNR = ABS( 6*Fc_t*dx1_t*dy1_t*dz1_t - sc_tnew**3 + (sc_tnew - dz1_t)**3 + &
                        (sc_tnew - dy1_t)**3 + (sc_tnew - dx1_t)**3)
        end if

        resNR = max(resNR,ABS(sc_tnew - sc_t))

        sc_t = sc_tnew
        itsNR = itsNR + 1

        end do

        write(1,*) 'time=',time,'its NR = ',itsNR, 'res NR = ', resNR

        end subroutine NewtonR

        ! Returns the volume of fluid in a given region (sub_x, sub_y, sub_z) within a cell
        ! The region must contain origin of s 
        function del_Vf(nor_xt, nor_yt, nor_zt, st, sub_x, sub_y, sub_z,phi_ijk)

        implicit none
        real :: nor_xt, nor_yt, nor_zt, st, sub_x, sub_y, sub_z, del_Vf
        real :: t1, t2, t3, sm, sc, Fc, dx1t, dy1t, dz1t, phi_ijk

        t1 = sub_x*ABS(nor_xt)
        t2 = sub_y*ABS(nor_yt)
        t3 = sub_z*ABS(nor_zt)

        ! Ensuring that dx1t >= dy1t >= dz1t        
        dx1t = max(t1,t2,t3)
        dz1t = min(t1,t2,t3)
        dy1t = t1 + t2 + t3 - dx1t - dz1t

        sm = dx1t + dy1t + dz1t
        sc = min(st,sm-st) ! st = s(i,j,k)

        ! For non-interfacial cells
        if (((nor_xt == 0.0) .and. (nor_yt == 0.0) .and. (nor_zt == 0.0)) .or. (ABS(phi_ijk) > eps)) then

                if (phi_ijk > 0) then
                        Fc = 1
                        del_Vf = Fc*sub_x*sub_y*sub_z
                else
                        Fc = 0
                        del_Vf = Fc*sub_x*sub_y*sub_z
                end if

        else
        ! For interfacial cells

        if ( sm <= st ) then
                Fc = 1
                del_Vf = Fc*sub_x*sub_y*sub_z
        else if ( st < 0 ) then
                Fc = 0
                del_Vf = Fc*sub_x*sub_y*sub_z
        else 
                if ( sc < dz1t ) then

                        Fc = (sc**3)/(6*dx1t*dy1t*dz1t) 

                else if ( ( sc >= dz1t ) .and. ( sc < dy1t ) ) then 

                        Fc = (sc**2 - dz1t*sc + ((dz1t**2)/3))/(2*dx1t*dy1t)

                else if ( ( sc >= (dy1t + dz1t) ) .and. ( sc <= dx1t ) ) then

                        Fc = (2*sc - dy1t - dz1t)/(2*dx1t) 

                else 

                        if (sc < dx1t) then
                                Fc = (sc**3 - (sc - dz1t)**3 - (sc - dy1t)**3)/(6*dx1t*dy1t*dz1t) 
                        else
                                Fc = (sc**3 - (sc - dz1t)**3 - (sc - dy1t)**3 - (sc - dx1t)**3)/(6*dx1t*dy1t*dz1t)
                        end if

                end if

        
                if (st <= 0.5*sm) then
                      del_Vf = Fc*sub_x*sub_y*sub_z  
                else
                      del_Vf = (1-Fc)*sub_x*sub_y*sub_z
                end if

        end if

        end if

        end function del_Vf
        


        ! To calculate F on i-1/2 and i+1/2 faces
        ! F_sweep = F (x-y-z sweep)/Fstar (z-x-y zweep)/Fstar2 (y-z-x sweep)
        subroutine FfaceCenter_xSweep(F_sweep,phi_sweep)

        implicit none
        integer :: i,j,k
        real, dimension(Nx+2, Ny+2, Nz+2),intent(in) :: F_sweep, phi_sweep
        real :: delta_Vf, ut, nxt, nyt, nzt, s_t
        ! ut = u(i,j,k), nxt = nor_x(i,j,k), nyt = nor_y(i,j,k), nzt = nor_z(i,j,k), s_t = s(i,j,k)        

        do k = 2,Nz+1,1
           do j = 2,Ny+1,1
              do i = 1,Nx+1,1

                ut = u(i,j,k)

                ! right face
                if (ut > 0) then
        
                   nxt = nor_x(i,j,k)
                   nyt = nor_y(i,j,k)
                   nzt = nor_z(i,j,k)
                   s_t = s(i,j,k)
                
                   if ( nxt*ut >= 0) then
        delta_Vf = del_Vf(nxt,nyt,nzt,s_t,del_t*ABS(ut),dye(j),dze(k), phi_sweep(i,j,k))                                        
                   else 
        delta_Vf = F_sweep(i,j,k)*dxe(i)*dye(j)*dze(k) - &
                   del_Vf(nxt,nyt,nzt,s_t,dxe(i) - del_t*ABS(ut),dye(j),dze(k),phi_sweep(i,j,k))
                   end if                
                        
                else if (ut < 0) then

                   nxt = nor_x(i+1,j,k)
                   nyt = nor_y(i+1,j,k)
                   nzt = nor_z(i+1,j,k)
                   s_t = s(i+1,j,k)

                   if ( nxt*ut >= 0) then
           delta_Vf = del_Vf(nxt,nyt,nzt,s_t,del_t*ABS(ut),dye(j),dze(k),phi_sweep(i+1,j,k))                                  
                   else
           delta_Vf = F_sweep(i+1,j,k)*dxe(i+1)*dye(j)*dze(k) - &
                      del_Vf(nxt,nyt,nzt,s_t,dxe(i+1) - del_t*ABS(ut),dye(j),dze(k),phi_sweep(i+1,j,k))   
                   end if
                                
                else                
                end if                                

                if (ut == 0) then
                        if (i == 1) then
                          F_ur(i,j,k) = F_sweep(i,j,k) + (dxe(i)/2)*((F_sweep(i+1,j,k) - &
                                          F_sweep(i,j,k))/(xc(i+1) - xc(i)))
                        else
                          F_ur(i,j,k) = F_sweep(i,j,k) + (dxe(i)/2)*((F_sweep(i+1,j,k) - &
                                          F_sweep(i-1,j,k))/(xc(i+1) - xc(i-1)))
                        end if

                else
                        F_ur(i,j,k) = delta_Vf/(ABS(ut)*del_t*dye(j)*dze(k))
                end if

                ! left face for next cell
                F_ul(i+1,j,k) = F_ur(i,j,k)

              end do
           end do
        end do

        end subroutine FfaceCenter_xSweep


        ! To calculate F on j-1/2 and j+1/2 faces
        ! F_sweep = F (y-z-x sweep)/Fstar (x-y-z zweep)/Fstar2 (z-x-y sweep)
        subroutine FfaceCenter_ySweep(F_sweep,phi_sweep)

        implicit none
        integer :: i,j,k
        real, dimension(Nx+2, Ny+2, Nz+2),intent(in) :: F_sweep,phi_sweep
        real :: delta_Vf, vt, nxt, nyt, nzt, s_t
        ! vt = v(i,j,k), nxt = nor_x(i,j,k), nyt = nor_y(i,j,k), nzt = nor_z(i,j,k), s_t = s(i,j,k)

        do k = 2,Nz+1,1
           do j = 1,Ny+1,1
              do i = 2,Nx+1,1

                vt = v(i,j,k)

                ! top face
                if (vt > 0) then

                   nxt = nor_x(i,j,k)
                   nyt = nor_y(i,j,k)
                   nzt = nor_z(i,j,k)
                   s_t = s(i,j,k)

                   if ( nyt*vt >= 0) then
                     delta_Vf = del_Vf(nxt,nyt,nzt,s_t,dxe(i),del_t*ABS(vt),dze(k),phi_sweep(i,j,k))
                   else
                     delta_Vf = F_sweep(i,j,k)*dxe(i)*dye(j)*dze(k) - &
                                del_Vf(nxt,nyt,nzt,s_t,dxe(i),dye(j) - del_t*ABS(vt),dze(k),phi_sweep(i,j,k))
                   end if

                else if (vt < 0) then

                   nxt = nor_x(i,j+1,k)
                   nyt = nor_y(i,j+1,k)
                   nzt = nor_z(i,j+1,k)
                   s_t = s(i,j+1,k)

                   if ( nyt*vt >= 0) then
                     delta_Vf = del_Vf(nxt,nyt,nzt,s_t,dxe(i),del_t*ABS(vt),dze(k),phi_sweep(i,j+1,k))
                   else
                     delta_Vf = F_sweep(i,j+1,k)*dxe(i)*dye(j+1)*dze(k) - &
                                del_Vf(nxt,nyt,nzt,s_t,dxe(i),dye(j+1) - del_t*ABS(vt),dze(k),phi_sweep(i,j+1,k))
                   end if

                else
                end if

                if (vt == 0) then
                        if (j == 1) then
                          F_vt(i,j,k) = F_sweep(i,j,k) + (dye(j)/2)*((F_sweep(i,j+1,k) - &
                                          F_sweep(i,j,k))/(yc(j+1) - yc(j)))
                        else
                          F_vt(i,j,k) = F_sweep(i,j,k) + (dye(j)/2)*((F_sweep(i,j+1,k) - &
                                          F_sweep(i,j-1,k))/(yc(j+1) - yc(j-1)))
                        end if

                else
                        F_vt(i,j,k) = delta_Vf/(ABS(vt)*del_t*dxe(i)*dze(k))
                end if

                ! bottom face for the next cell
                F_vb(i,j+1,k) = F_vt(i,j,k)

              end do
           end do
        end do

        end subroutine FfaceCenter_ySweep


        ! To calculate F on k-1/2 and k+1/2 faces
        ! F_sweep = F (z-x-y sweep)/Fstar (y-z-x zweep)/Fstar2 (x-y-z sweep)
        subroutine FfaceCenter_zSweep(F_sweep,phi_sweep)

        implicit none
        integer :: i,j,k
        real, dimension(Nx+2, Ny+2, Nz+2),intent(in) :: F_sweep,phi_sweep
        real :: delta_Vf, wt, nxt, nyt, nzt, s_t
        ! wt = w(i,j,k), nxt = nor_x(i,j,k), nyt = nor_y(i,j,k), nzt =  nor_z(i,j,k), s_t = s(i,j,k)

        do k = 1,Nz+1,1
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1

                wt = w(i,j,k)

                ! south face
                if (wt > 0) then

                   nxt = nor_x(i,j,k)
                   nyt = nor_y(i,j,k)
                   nzt = nor_z(i,j,k)
                   s_t = s(i,j,k)

                   if ( nzt*wt >= 0) then
                     delta_Vf = del_Vf(nxt,nyt,nzt,s_t,dxe(i),dye(j),del_t*ABS(wt),phi_sweep(i,j,k))
                   else
                     delta_Vf = F_sweep(i,j,k)*dxe(i)*dye(j)*dze(k) - &
                                del_Vf(nxt,nyt,nzt,s_t,dxe(i),dye(j),dze(k) - del_t*ABS(wt),phi_sweep(i,j,k))
                   end if

                else if (wt < 0) then

                   nxt = nor_x(i,j,k+1)
                   nyt = nor_y(i,j,k+1)
                   nzt = nor_z(i,j,k+1)
                   s_t = s(i,j,k+1)

                   if ( nzt*wt >= 0) then
                     delta_Vf = del_Vf(nxt,nyt,nzt,s_t,dxe(i),dye(j),del_t*ABS(wt),phi_sweep(i,j,k+1))
                   else
                     delta_Vf = F_sweep(i,j,k+1)*dxe(i)*dye(j)*dze(k+1) - &
                                del_Vf(nxt,nyt,nzt,s_t,dxe(i),dye(j),dze(k+1) - del_t*ABS(wt),phi_sweep(i,j,k+1))
                   end if

                else
                end if

                if (wt == 0) then
                        if (k == 1) then
                          F_ws(i,j,k) = F_sweep(i,j,k) + (dze(k)/2)*((F_sweep(i,j,k+1) - &
                                          F_sweep(i,j,k))/(zc(k+1) - zc(k)))
                        else
                          F_ws(i,j,k) = F_sweep(i,j,k) + (dze(k)/2)*((F_sweep(i,j,k+1) - &
                                          F_sweep(i,j,k-1))/(zc(k+1) - zc(k-1)))
                        end if

                else
                        F_ws(i,j,k) = delta_Vf/(ABS(wt)*del_t*dxe(i)*dye(j))
                end if

                ! north face for the next cell
                F_wn(i,j,k+1) = F_ws(i,j,k)

              end do
           end do
        end do

        end subroutine FfaceCenter_zSweep


        ! Reinitialisation of LS function 
        ! F_tem = F_np1, phi_tem = phi_np1
        ! Corrects phi_tem according to the reconstructed interface
        subroutine reinitialisePhi(F_tem, phi_tem)

        real, dimension(Nx+2, Ny+2, Nz+2) :: d_cor, phi_tt
        real, dimension(Nx+2, Ny+2, Nz+2), intent(in) :: F_tem
        real, dimension(Nx+2, Ny+2, Nz+2), intent(out) :: phi_tem
        integer :: i,j,k,Kk,i_prime,j_prime,k_prime, i_p,j_p,k_p,l,m,n
        integer :: inew, jnew,knew, iin, jjn, kkn
        real, dimension(3) :: xv,xp,x_prime, ps1, ps2, pis, pr1, pr2, tmp1,tmp2
        real :: d_prod, nnx, nny, nnz, d_max, nxt, nyt, nzt, nxc, nyc, nzc, dt, dc, nor, var1, var2
       
        phi_aux  = 10000 
        phi_tt = 0
        min_d = 10000
        ! phi is determined accurately in (i-Kk,j-Kk,k-Kk) to (i+Kk,j+Kk,k+Kk) region surrounding interfacial cell (i,j,k)
        Kk = 4
        d_cor = 0
        d_prime = 0


        tmp1(1) = xe(1)
        tmp1(2) = ye(1)
        tmp1(3) = ze(1)
        tmp2(1) = xe(2)
        tmp2(2) = ye(2)
        tmp2(3) = ze(2)

        ! Maximum diagonal length of all computational cells
        d_max = Euc(tmp1,tmp2)
        
        ! d_cor is the corrected d obtained from the corrected s
        do k = 2,Nz+1,1
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1
                d_cor(i,j,k) = s(i,j,k) - 0.5*( ABS(nor_x(i,j,k))*dxe(i) + &
                               ABS(nor_y(i,j,k))*dye(j) + ABS(nor_z(i,j,k))*dze(k) )
              end do
           end do
        end do

        ! d_prime is the distance between the reconstructed interface and origin (0,0,0)
        do k = 2,Nz+1,1
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1
                d_prime(i,j,k) = d_cor(i,j,k) - ( nor_x(i,j,k)*xc(i) + &
                                 nor_y(i,j,k)*yc(j) + nor_z(i,j,k)*zc(k) )
              end do
           end do
        end do

        do k = 2,Nz+1,1
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1

		flag = 0

                 ! Finding cells where interface exists
                 if ( (F_tem(i,j,k) > F_low) .and. (F_tem(i,j,k) < F_up) .and. (ABS(phi_tem(i,j,k)) <= eps) ) then
                    do k_prime = k-1,k+1,1
                      do j_prime = j-1,j+1,1
                        do i_prime = i-1,i+1,1
                           if ( ( (phi_tem(i,j,k)*phi_tem(i_prime,j_prime,k_prime) < 0) &
                              .or. (phi_tem(i,j,k) == 0) ) .and. ((i-i_prime)*(j-j_prime)*(k-k_prime) <= 1) ) then

				flag = 1

                                nnx = nor_x(i,j,k)
                                nny = nor_y(i,j,k)
                                nnz = nor_z(i,j,k)

                                do k_p = k-Kk,k+Kk,1
                                  do j_p = j-Kk,j+Kk,1
                                    do i_p = i-Kk,i+Kk,1

                                  ! Checking that the cells near the interfacial cell (within a band of Kk) indeed are inside the
                                  ! computational domain  
                                  if ( (i_p >= 1) .and. (i_p <= Nx+2) .and. (j_p >= 1) .and. &
                                       (j_p <= Ny+2) .and. (k_p >= 1) .and. (k_p <= Nz+2) ) then  
                                        
                                        ! x_prime = cell center (i_p,j_p,k_p) vector
                                        x_prime(1) = xc(i_p)
                                        x_prime(2) = yc(j_p)
                                        x_prime(3) = zc(k_p)

                                        ! when i_p = i, j_p = j, k_p = k (interfacial cell)
                                        if ((i_p == i) .and. (j_p == j) .and. (k_p == k)) then
                                           phi_aux(i_p,j_p,k_p) = SIGN(one,F_tem(i_p,j_p,k_p)-0.5)*ABS(Dis(i,j,k,x_prime))
                                        else
                                          
                                           ! Determining xv = point on the  boundary of cell (i,j,k) with the shortest distance to cell center (i_p,j_p,k_p)
                                           l = max(-1, min(1, i_p - i))
                                           m = max(-1, min(1, j_p - j))
                                           n = max(-1, min(1, k_p - k))

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

                                           if (n == -1) then
                                                xv(3) = ze(k-1)
                                           else if (n == 0) then
                                                xv(3) = zc(k)
                                           else if (n == 1) then
                                                xv(3) = ze(k)
                                           else
                                           end if

                                           if ( Dis(i,j,k,xv)*SIGN(one,F_tem(i_p,j_p,k_p)-0.5) <= 0 ) then
                                              phi_aux(i_p,j_p,k_p) = &
                                              SIGN(one,F_tem(i_p,j_p,k_p)-0.5)*min( Euc(xv,x_prime),ABS(phi_aux(i_p,j_p,k_p)) )
                                           else
                                           end if
                                                
                                              ! xp = projection of x_prime onto the interface in cell (i,j,k)
                                              xp(1) = x_prime(1) - Dis(i,j,k,x_prime)*nnx  
                                              xp(2) = x_prime(2) - Dis(i,j,k,x_prime)*nny
                                              xp(3) = x_prime(3) - Dis(i,j,k,x_prime)*nnz
                                              ! Checking if xp is inside cell (i,j,k)
                                              if (((xp(1) > xe(i-1)) .and. (xp(1) < xe(i))) .and. &
                                                 ( (xp(2) > ye(j-1)) .and. (xp(2) < ye(j))) .and. &
                                                 ( (xp(3) > ze(k-1)) .and. (xp(3) < ze(k))) ) then

                                                phi_aux(i_p,j_p,k_p) = &
                                                SIGN(one,F_tem(i_p,j_p,k_p)-0.5)*min( Euc(xp,x_prime),ABS(phi_aux(i_p,j_p,k_p)) )

                                              else
                                    
                                                ! (Step 4)
                                                min_d = 10000

                                                do knew = k-1,k,1
                                                 do jnew = j-1,j,1
                                                  do inew = i-1,i,1

                                                    do kkn = k-1,k,1
                                                     do jjn = j-1,j,1
                                                      do iin = i-1,i,1

                                                        ps1(1) = xe(inew)
                                                        ps1(2) = ye(jnew)
                                                        ps1(3) = ze(knew)
                                                        ps2(1) = xe(iin)
                                                        ps2(2) = ye(jjn)
                                                        ps2(3) = ze(kkn)

                                                        d_prod = Dis(i,j,k,ps1)*Dis(i,j,k,ps2)
                                                        
                                                        ! edge corners separated by dze(i)        
                                                        if ( (knew < kkn) .and. (jnew == jjn) .and. (inew == iin)) then
                                                          if (d_prod > 0) then
                                                          else
                                                            pis(1) = ps1(1) 
                                                            pis(2) = ps1(2)
                                                            pis(3) = (nnx*pis(1) + nny*pis(2) + d_prime(i,j,k))/(-nnz)
                                                            min_d = min( Euc(pis,x_prime),min_d ) 
   
                                                            nor = SQRT(nnx**2 + nnz**2)
                                                            nyt = nny/nor
                                                            nzt = nnz/nor
                                                            dt  = (d_prime(i,j,k) + nnx*pis(1))/nor
                                                            var1 = nyt*x_prime(2) + nzt*x_prime(3) + dt
                                                            nor = SQRT(nny**2 + nnz**2)
                                                            nxc = nnx/nor
                                                            nzc = nnz/nor
                                                            dc  = (d_prime(i,j,k) + nny*pis(2))/nor
                                                            var2 = nxc*x_prime(1) + nzc*x_prime(3) + dc

                                                            pr1(1) = pis(1)
                                                            pr1(2) = x_prime(2) - var1*nyt
                                                            pr1(3) = x_prime(3) - var1*nzt

                                                            pr2(1) = x_prime(1) - var2*nxc
                                                            pr2(2) = pis(2)
                                                            pr2(3) = x_prime(3) - var2*nzc 

                                                            ! Checking if pr1 is inside cell (i,j,k)
                                              if (((pr1(1) >= xe(i-1)) .and. (pr1(1) <= xe(i))) .and. &
                                                 ( (pr1(2) >= ye(j-1)) .and. (pr1(2) <= ye(j))) .and. &
                                                 ( (pr1(3) >= ze(k-1)) .and. (pr1(3) <= ze(k))) ) then

                                                min_d = &
                                                min( Euc(pr1,x_prime), min_d )

                                              else
                                              end if
                                                              ! Checking if pr2 is inside cell (i,j,k)
                                              if (((pr2(1) >= xe(i-1)) .and. (pr2(1) <= xe(i))) .and. &
                                                 ( (pr2(2) >= ye(j-1)) .and. (pr2(2) <= ye(j))) .and. &
                                                 ( (pr2(3) >= ze(k-1)) .and. (pr2(3) <= ze(k))) ) then

                                                min_d = &
                                                min( Euc(pr2,x_prime), min_d )

                                              else
                                              end if

                                                            end if     

                                                        ! edge corners separated by dye(j)         
                                                        else if ( (knew == kkn) .and. (jnew < jjn) .and. (inew == iin)) then
                                                          if (d_prod > 0) then
                                                          else
                                                            pis(1) = ps1(1)
                                                            pis(3) = ps1(3)    
                                                            pis(2) = (nnx*pis(1) + nnz*pis(3) + d_prime(i,j,k))/(-nny)

                                                            min_d = min( Euc(pis,x_prime),min_d )

                                                            nor = SQRT(nnx**2 + nny**2)
                                                            nyt = nny/nor
                                                            nzt = nnz/nor
                                                            dt  = (d_prime(i,j,k) + nnx*pis(1))/nor
                                                            var1 = nyt*x_prime(2) + nzt*x_prime(3) + dt
                                                            nor = SQRT(nny**2 + nnz**2)
                                                            nxc = nnx/nor
                                                            nyc = nny/nor
                                                            dc  = (d_prime(i,j,k) + nnz*pis(3))/nor
                                                            var2 = nxc*x_prime(1) + nyc*x_prime(2) + dc

                                                            pr1(1) = pis(1)
                                                            pr1(2) = x_prime(2) - var1*nyt
                                                            pr1(3) = x_prime(3) - var1*nzt

                                                            pr2(1) = x_prime(1) - var2*nxc
                                                            pr2(2) = x_prime(2) - var2*nyc 
                                                            pr2(3) = pis(3)

                                                            ! Checking if pr1 is
                                                            ! inside cell
                                                            ! (i,j,k)
                                              if (((pr1(1) >= xe(i-1)) .and. (pr1(1) <= xe(i))) .and. &
                                                 ( (pr1(2) >= ye(j-1)) .and. (pr1(2) <= ye(j))) .and. &
                                                 ( (pr1(3) >= ze(k-1)) .and. (pr1(3) <= ze(k))) ) then

                                                min_d = &
                                                min( Euc(pr1,x_prime), ABS(min_d) )

                                              else
                                              end if
                                                              ! Checking if pr2
                                                              ! is inside cell
                                                              ! (i,j,k)
                                              if (((pr2(1) >= xe(i-1)) .and. (pr2(1) <= xe(i))) .and. &
                                                 ( (pr2(2) >= ye(j-1)) .and. (pr2(2) <= ye(j))) .and. &
                                                 ( (pr2(3) >= ze(k-1)) .and. (pr2(3) <= ze(k))) ) then

                                                min_d = &
                                                min( Euc(pr2,x_prime), ABS(min_d) )

                                              else
                                              end if

                                                          end if
                                                         
                                                        ! edge corners separated by dxe(i)
                                                        else if ( (knew == kkn) .and. (jnew == jjn) .and. (inew < iin)) then
                                                          if (d_prod > 0) then
                                                          else
                                                            pis(2) = ps1(2)
                                                            pis(3) = ps1(3)
                                                            pis(1) = (nny*pis(2) + nnz*pis(3) + d_prime(i,j,k))/(-nnx)

                                                            min_d = min( Euc(pis,x_prime),min_d )

                                                            nor = SQRT(nnx**2 + nny**2)
                                                            nxt = nny/nor
                                                            nzt = nnz/nor
                                                            dt  = (d_prime(i,j,k) + nny*pis(2))/nor
                                                            var1 = nxt*x_prime(1) + nzt*x_prime(3) + dt
                                                            nor = SQRT(nnx**2 + nnz**2)
                                                            nxc = nnx/nor
                                                            nyc = nny/nor
                                                            dc  = (d_prime(i,j,k) + nnz*pis(3))/nor
                                                            var2 = nxc*x_prime(1) + nyc*x_prime(2) + dc

                                                            pr1(1) = x_prime(1) - var1*nxt 
                                                            pr1(2) = pis(2)
                                                            pr1(3) = x_prime(3) - var1*nzt

                                                            pr2(1) = x_prime(1) - var2*nxc
                                                            pr2(2) = x_prime(2) - var2*nyc
                                                            pr2(3) = pis(3)

                                                            ! Checking if pr1 is
                                                            ! inside cell
                                                            ! (i,j,k)
                                              if (((pr1(1) >= xe(i-1)) .and. (pr1(1) <= xe(i))) .and. &
                                                 ( (pr1(2) >= ye(j-1)) .and. (pr1(2) <= ye(j))) .and. &
                                                 ( (pr1(3) >= ze(k-1)) .and. (pr1(3) <= ze(k))) ) then

                                                min_d = &
                                                min( Euc(pr1,x_prime), min_d )

                                              else
                                              end if
                                                              ! Checking if pr2
                                                              ! is inside cell
                                                              ! (i,j,k)
                                              if (((pr2(1) >= xe(i-1)) .and. (pr2(1) <= xe(i))) .and. &
                                                 ( (pr2(2) >= ye(j-1)) .and. (pr2(2) <= ye(j))) .and. &
                                                 ( (pr2(3) >= ze(k-1)) .and. (pr2(3) <= ze(k))) ) then

                                                min_d = &
                                                min( Euc(pr2,x_prime), min_d )

                                              else
                                              end if


                                                          end if
                                                        else
                                                        end if

                                                      end do
                                                     end do
                                                    end do
        
                                                  end do
                                                 end do
                                                end do
                                                
                                                ! All the 12 edges have been
                                                ! checked for intersection with
                                                ! the interfacial plane
                                                phi_aux(i_p,j_p,k_p) = &
                                                SIGN(one,F_tem(i_p,j_p,k_p)-0.5)*min( ABS(min_d),ABS(phi_aux(i_p,j_p,k_p)) )
                                                       
                                              end if 
                                                     
                                        end if       

                                else
                                end if

                                    end do
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

			 if (flag == 1) then
                         exit
                        else
                        end if

                    end do

                 else
                 end if

              end do
           end do
        end do


        do k = 1,Nz+2,1
          do j = 1,Ny+2,1
             do i = 1,Nx+2,1

                if (phi_tem(i,j,k) > 0) then
                        phi_tt(i,j,k) = (Kk+1)*d_max
                else if (phi_tem(i,j,k) < 0) then
                        phi_tt(i,j,k) = -(Kk+1)*d_max
                else
                end if

             end do
          end do
        end do


        do k = 2,Nz+1,1
           do j = 2,Ny+1,1
              do i = 2,Nx+1,1

		flag = 0

                 ! Finding cells where interface exists
                 if ( (F_tem(i,j,k) > F_low) .and. (F_tem(i,j,k) < F_up) .and. (ABS(phi_tem(i,j,k)) <= eps) ) then
                    do k_prime = k-1,k+1,1
                      do j_prime = j-1,j+1,1
                        do i_prime = i-1,i+1,1
                           if ( ( (phi_tem(i,j,k)*phi_tem(i_prime,j_prime,k_prime) < 0) &
                              .or. (phi_tem(i,j,k) == 0) ) .and. ((i-i_prime)*(j-j_prime)*(k-k_prime) <= 1) ) then

				flag = 1

                                do k_p = k-Kk,k+Kk,1
                                  do j_p = j-Kk,j+Kk,1
                                    do i_p = i-Kk,i+Kk,1

                                    ! Checking that the cells near the interfacial cell (within a band of Kk) indeed are inside the
                                    ! computational domain
                                    if ( (i_p >= 1) .and. (i_p <= Nx+2) .and. (j_p >= 1) .and. &
                                       (j_p <= Ny+2) .and. (k_p >= 1) .and. (k_p <= Nz+2) ) then

                                      phi_tt(i_p,j_p,k_p) = phi_aux(i_p,j_p,k_p)

                                    else
                                    end if

                                    end do
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

			 if (flag == 1) then
                         exit
                        else
                        end if

                    end do

                 else
                 end if

              end do
           end do
        end do

        phi_tem = phi_tt

        end subroutine reinitialisePhi

        ! Returns distance of point xd (vector) from the plane in cell (ii,jj,kk)
        function Dis(ii,jj,kk,xd)
                real :: Dis
                integer :: ii,jj,kk
                real, dimension(3) :: xd

                Dis = nor_x(ii,jj,kk)*xd(1) + nor_y(ii,jj,kk)*xd(2) + &
                      nor_z(ii,jj,kk)*xd(3) + d_prime(ii,jj,kk)
        end function Dis

        ! Returns Euclidean distance between points p1 and p2 (vectors)
        function Euc(p1,p2)
                real :: Euc
                real, dimension(3) :: p1,p2
                Euc = ( (p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2 )**0.5
        end function Euc


        ! To calculate zeta_star (= Fstar .or. phi_star)
        ! zeta_ul = zeta_i-1/2 or zeta_j-1/2 or zeta_k-1/2 (depending on the sweep direction)
        ! zeta_ur = zeta_i+1/2 or zeta_j+1/2 or zeta_k+1/2 (zeta = F or phi)
        subroutine calcStar(zeta, zeta_ul, zeta_ur, zeta_star, sweepd)

        implicit none
        integer :: i,j,k, sweepd
        real, dimension(Nx+2, Ny+2, Nz+2), intent(in) :: zeta, zeta_ul, zeta_ur
        real, dimension(Nx+2, Ny+2, Nz+2), intent(out) :: zeta_star
                
        ! sweepd = 1 for x-sweep, 2 for y-sweep and 3 for z-sweep

        if (sweepd == 1) then

                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_star(i,j,k) = ( zeta(i,j,k) - (del_t/dxe(i))*( zeta_ur(i,j,k)*u(i,j,k) - &
                                  zeta_ul(i,j,k)*u(i-1,j,k) ) )/( 1 - (del_t/dxe(i))*(u(i,j,k) - u(i-1,j,k)) ) 
                                end do
                        end do
                end do

        else if (sweepd == 2) then

                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_star(i,j,k) = ( zeta(i,j,k) - (del_t/dye(j))*( zeta_ur(i,j,k)*v(i,j,k) - &
                                  zeta_ul(i,j,k)*v(i,j-1,k) ) )/( 1 - (del_t/dye(j))*(v(i,j,k) - v(i,j-1,k)) )
                                end do
                        end do
                end do

        else 

                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_star(i,j,k) = ( zeta(i,j,k) - (del_t/dze(k))*( zeta_ur(i,j,k)*w(i,j,k) - &
                                  zeta_ul(i,j,k)*w(i,j,k-1) ) )/( 1 - (del_t/dze(k))*(w(i,j,k) - w(i,j,k-1)) )
                                end do
                        end do
                end do

        end if

        ! Boundary Conditions (Neumann)
        ! x = 0 (Inlet)
        zeta_star(1,:,:)    = zeta_star(2,:,:)

        ! x = Lx (Outlet)
        zeta_star(Nx+2,:,:) = zeta_star(Nx+1,:,:)

        ! y = 0 (ground)
        zeta_star(:,1,:)    = zeta_star(:,2,:)

        ! y = Ly (top)
        zeta_star(:,Ny+2,:) = zeta_star(:,Ny+1,:)

        ! z = 0 (front)
        zeta_star(:,:,1)    = zeta_star(:,:,2)

        ! z = Lz (back)
        zeta_star(:,:,Nz+2) = zeta_star(:,:,Nz+1)

        end subroutine calcStar

        ! To calculate zeta_star2 (= Fstar2 .or. phi_star2)
        ! zeta_vb = zeta_i-1/2 or zeta_j-1/2 or zeta_k-1/2 (depending on the sweep direction)
        ! zeta_vt = zeta_i+1/2 or zeta_j+1/2 or zeta_k+1/2 (zeta = F or phi)
        subroutine calcStar2(zeta_star, zeta_vb, zeta_vt, zeta_star2, sweepd)

        implicit none
        integer :: i,j,k, sweepd
        real, dimension(Nx+2, Ny+2, Nz+2), intent(in) :: zeta_star, zeta_vb, zeta_vt
        real, dimension(Nx+2, Ny+2, Nz+2), intent(out) :: zeta_star2

        ! sweepd = 1 for x-sweep, 2 for y-sweep and 3 for z-sweep

        if (sweepd == 1) then

                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_star2(i,j,k) =  zeta_star(i,j,k)*( 1 + (del_t/dxe(i))*(u(i,j,k) - u(i-1,j,k)) ) - &
                                  (del_t/dxe(i))*( zeta_vt(i,j,k)*u(i,j,k) - zeta_vb(i,j,k)*u(i-1,j,k) )
                                end do
                        end do
                end do


        else if (sweepd == 2) then

                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_star2(i,j,k) =  zeta_star(i,j,k)*( 1 + (del_t/dye(j))*(v(i,j,k) - v(i,j-1,k)) ) - &
                                  (del_t/dye(j))*( zeta_vt(i,j,k)*v(i,j,k) - zeta_vb(i,j,k)*v(i,j-1,k) ) 
                                end do
                        end do
                end do

        else

                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_star2(i,j,k) =  zeta_star(i,j,k)*( 1 + (del_t/dze(k))*(w(i,j,k) - w(i,j,k-1)) ) - &
                                  (del_t/dze(k))*( zeta_vt(i,j,k)*w(i,j,k) - zeta_vb(i,j,k)*w(i,j,k-1) )
                                end do
                        end do
                end do


        end if

        ! Boundary Conditions (Neumann)
        ! x = 0 (Inlet)
        zeta_star2(1,:,:)    = zeta_star2(2,:,:)

        ! x = Lx (Outlet)
        zeta_star2(Nx+2,:,:) = zeta_star2(Nx+1,:,:)

        ! y = 0 (ground)
        zeta_star2(:,1,:)    = zeta_star2(:,2,:)

        ! y = Ly (top)
        zeta_star2(:,Ny+2,:) = zeta_star2(:,Ny+1,:)

        ! z = 0 (front)
        zeta_star2(:,:,1)    = zeta_star2(:,:,2)

        ! z = Lz (back)
        zeta_star2(:,:,Nz+2) = zeta_star2(:,:,Nz+1)

        end subroutine calcStar2

        ! To calculate zeta_np1 (= F_np1 .or. phi_np1) 
        ! zeta_wn = zeta_i-1/2 or zeta_j-1/2 or zeta_k-1/2 (depending on the sweep direction)
        ! zeta_ws = zeta_i+1/2 or zeta_j+1/2 or zeta_k+1/2 (zeta = F or phi)
        subroutine calc_np1(zeta_star, zeta_star2, zeta_wn, zeta_ws, zeta_np1, sweepd)

        implicit none
        integer :: i,j,k, sweepd
        real, dimension(Nx+2, Ny+2, Nz+2), intent(in) :: zeta_star, zeta_star2, zeta_wn, zeta_ws
        real, dimension(Nx+2, Ny+2, Nz+2), intent(out) :: zeta_np1

        ! sweepd = 1 for x-sweep, 2 for y-sweep and 3 for z-sweep

        if (sweepd == 1) then
        
                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_np1(i,j,k) =  zeta_star2(i,j,k) + zeta_star(i,j,k)*(del_t/dxe(i))*(u(i,j,k) - &
                                  u(i-1,j,k)) - (del_t/dxe(i))*( zeta_ws(i,j,k)*u(i,j,k) - zeta_wn(i,j,k)*u(i-1,j,k) )
                                end do
                        end do
                end do

        else if (sweepd == 2) then

                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_np1(i,j,k) =  zeta_star2(i,j,k) + zeta_star(i,j,k)*(del_t/dye(j))*(v(i,j,k) - &
                                  v(i,j-1,k)) - (del_t/dye(j))*( zeta_ws(i,j,k)*v(i,j,k) - zeta_wn(i,j,k)*v(i,j-1,k) )
                                end do
                        end do
                end do

        else
           
                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                  zeta_np1(i,j,k) =  zeta_star2(i,j,k) + zeta_star(i,j,k)*(del_t/dze(k))*(w(i,j,k) - &
                                  w(i,j,k-1)) - (del_t/dze(k))*( zeta_ws(i,j,k)*w(i,j,k) - zeta_wn(i,j,k)*w(i,j,k-1) )
                                end do
                        end do
                end do
        
        end if

        ! Boundary Conditions (Neumann)
        ! x = 0 (Inlet)
        zeta_np1(1,:,:)    = zeta_np1(2,:,:)

        ! x = Lx (Outlet)
        zeta_np1(Nx+2,:,:) = zeta_np1(Nx+1,:,:)

        ! y = 0 (ground)
        zeta_np1(:,1,:)    = zeta_np1(:,2,:)

        ! y = Ly (top)
        zeta_np1(:,Ny+2,:) = zeta_np1(:,Ny+1,:)

        ! z = 0 (front)
        zeta_np1(:,:,1)    = zeta_np1(:,:,2)

        ! z = Lz (back)
        zeta_np1(:,:,Nz+2) = zeta_np1(:,:,Nz+1)

        end subroutine calc_np1

        ! To calculate phi on i-1/2 and i+1/2 faces
        ! phi_sweep = phi (x-y-z sweep)/phi_star (z-x-y zweep)/phi_star2 (y-z-x sweep)
        subroutine phiFaceCenter_xSweep(phi_sweep)

        implicit none
        integer :: i,j,k
        real, dimension(Nx+2, Ny+2, Nz+2) :: phi_sweep

        do k = 2,Nz+1,1
                do j = 2,Ny+1,1
                        do i = 1,Nx+1,1

                                ! right face
                                if (u(i,j,k) >= 0) then

                                        if (i == 1) then
                                                phi_ur(i,j,k) = phi_sweep(i,j,k) + (dxe(i)/2)*(1 - &
                                                u(i,j,k)*(del_t/dxe(i)))*((phi_sweep(i+1,j,k) - &
                                                phi_sweep(i,j,k))/(xc(i+1) - xc(i)))
                                        else
                                                phi_ur(i,j,k) = phi_sweep(i,j,k) + (dxe(i)/2)*(1 - &
                                                u(i,j,k)*(del_t/dxe(i)))*((phi_sweep(i+1,j,k) - &
                                                phi_sweep(i-1,j,k))/(xc(i+1) - xc(i-1)))
                                        end if

                                else 

                                        if (i == Nx+1) then
                                                phi_ur(i,j,k) = phi_sweep(i,j,k) + (dxe(i)/2)*(1 + &
                                                u(i,j,k)*(del_t/dxe(i)))*((phi_sweep(i+1,j,k) - &
                                                phi_sweep(i,j,k))/(xc(i+1) - xc(i)))
                                        else
                                                phi_ur(i,j,k) = phi_sweep(i,j,k) + (dxe(i)/2)*(1 + &
                                                u(i,j,k)*(del_t/dxe(i)))*((phi_sweep(i+2,j,k) - &
                                                phi_sweep(i,j,k))/(xc(i+2) - xc(i)))
                                        end if

                                end if
                        
                                ! left face
                                phi_ul(i+1,j,k) = phi_ur(i,j,k)                        
                        end do
                end do
        end do

        end subroutine phiFaceCenter_xSweep

        ! To calculate phi on j-1/2 and j+1/2 faces
        ! phi_sweep = phi (y-z-x sweep)/phi_star (x-y-z zweep)/phi_star2 (z-x-y  sweep)
        subroutine phiFaceCenter_ySweep(phi_sweep)

        implicit none
        integer :: i,j,k
        real, dimension(Nx+2, Ny+2, Nz+2) :: phi_sweep

        do k = 2,Nz+1,1
                do j = 1,Ny+1,1
                        do i = 2,Nx+1,1

                                ! top face
                                if (v(i,j,k) >= 0) then

                                        if (j == 1) then
                                                phi_vt(i,j,k) = phi_sweep(i,j,k) + (dye(j)/2)*(1 - &
                                                v(i,j,k)*(del_t/dye(j)))*((phi_sweep(i,j+1,k) - &
                                                phi_sweep(i,j,k))/(yc(j+1) - yc(j)))
                                        else
                                                phi_vt(i,j,k) = phi_sweep(i,j,k) + (dye(j)/2)*(1 - &
                                                v(i,j,k)*(del_t/dye(j)))*((phi_sweep(i,j+1,k) - &
                                                phi_sweep(i,j-1,k))/(yc(j+1) - yc(j-1)))
                                        end if

                                else 

                                        if (j == Ny+1) then
                                                phi_vt(i,j,k) = phi_sweep(i,j,k) + (dye(j)/2)*(1 + &
                                                v(i,j,k)*(del_t/dye(j)))*((phi_sweep(i,j+1,k) - &
                                                phi_sweep(i,j,k))/(yc(j+1) - yc(j)))
                                        else
                                                phi_vt(i,j,k) = phi_sweep(i,j,k) + (dye(j)/2)*(1 + &
                                                v(i,j,k)*(del_t/dye(j)))*((phi_sweep(i,j+2,k) - &
                                                phi_sweep(i,j,k))/(yc(j+2) - yc(j)))
                                        end if

                                end if

                                ! bottom face
                                phi_vb(i,j+1,k) = phi_vt(i,j,k)

                        end do
                end do
        end do

        end subroutine phiFaceCenter_ySweep

        ! To calculate phi on k-1/2 and k+1/2 faces
        ! phi_sweep = phi (z-x-y sweep)/phi_star (y-z-x zweep)/phi_star2 (x-y-z sweep)
        subroutine phiFaceCenter_zSweep(phi_sweep)

        implicit none
        integer :: i,j,k
        real, dimension(Nx+2, Ny+2, Nz+2) :: phi_sweep

        do k = 1,Nz+1,1
                do j = 2,Ny+1,1
                        do i = 2,Nx+1,1

                                ! south face
                                if (w(i,j,k) >= 0) then

                                        if (k == 1) then
                                                phi_ws(i,j,k) = phi_sweep(i,j,k) + (dze(k)/2)*(1 - &
                                                w(i,j,k)*(del_t/dze(k)))*((phi_sweep(i,j,k+1) - &
                                                phi_sweep(i,j,k))/(zc(k+1) - zc(k)))
                                        else
                                                phi_ws(i,j,k) = phi_sweep(i,j,k) + (dze(k)/2)*(1 - &
                                                w(i,j,k)*(del_t/dze(k)))*((phi_sweep(i,j,k+1) - &
                                                phi_sweep(i,j,k-1))/(zc(k+1) - zc(k-1)))
                                        end if

                                else 

                                        if (k == Nz+1) then
                                                phi_ws(i,j,k) = phi_sweep(i,j,k) + (dze(k)/2)*(1 + &
                                                w(i,j,k)*(del_t/dze(k)))*((phi_sweep(i,j,k+1) - &
                                                phi_sweep(i,j,k))/(zc(k+1) - zc(k)))
                                        else
                                                phi_ws(i,j,k) = phi_sweep(i,j,k) + (dze(k)/2)*(1 + &
                                                w(i,j,k)*(del_t/dze(k)))*((phi_sweep(i,j,k+2) - &
                                                phi_sweep(i,j,k))/(zc(k+2) - zc(k)))
                                        end if

                                end if

                                ! north face
                                phi_wn(i,j,k+1) = phi_ws(i,j,k)

                        end do
                end do
        end do

        end subroutine phiFaceCenter_zSweep

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
                integer :: ikea, jkea, kkea
                real :: delt, delt_x, delt_y, delt_z, delx, dely, delz, umax, vmax, wmax
                umax = 0
                vmax = 0
                wmax = 0

                delx = minval(dxe(2:Nx+2))
                dely = minval(dye(2:Ny+2))
                delz = minval(dze(2:Nz+2))

                do kkea = 2,Nz+1,1
                  do jkea = 2,Ny+1,1
                    do ikea = 2,Nx+1,1
                        
                        if (umax < ABS(u(ikea,jkea,kkea))) then
                                umax = ABS(u(ikea,jkea,kkea))
                        else
                        end if

                        if (vmax < ABS(v(ikea,jkea,kkea))) then
                                vmax = ABS(v(ikea,jkea,kkea))
                        else
                        end if

                        if (wmax < ABS(w(ikea,jkea,kkea))) then
                                wmax = ABS(w(ikea,jkea,kkea))
                        else
                        end if

                    end do
                  end do
                end do

                delt_x = (0.25*delx)/umax
                delt_y = (0.25*dely)/vmax
                delt_z = (0.25*delz)/wmax

                if ((umax < (0.25*delx)) .and. (vmax < 0.25*dely) .and. (wmax < 0.25*delz)) then
                        delt = 0.01
                else        
                        delt = min(delt_x, delt_y, delt_z)
                end if

        end subroutine delta_t



end program main
