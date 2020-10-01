! ME634A Assignment-4 : 3D lid-driven cavity flow using RK3-CN method
! Code by Shreya Agrawal (160662) in Fortran 95

program main

     implicit none

     real :: Lx, Ly, Lz, Re, Lo, Uo, nu, del_t, resd, t_resdu, t_resdv, t_resdw, t_resdp, t_tol, time         
     real, dimension(:,:,:), allocatable :: u, v, w, p, u_np1, v_np1, w_np1, p_np1    ! velocity and pressure at present time and
                                                                                      ! the next time-step
     real, dimension(:,:,:), allocatable :: ustar, vstar, wstar, u1, v1, w1, p1, ustar2, vstar2, wstar2, u2, v2, w2, p2, &
                                            ustar3, vstar3, wstar3, qu0, qv0, qw0, qu1, qv1, qw1, qu2, qv2, qw2, qu3, qv3, qw3, &
                                            rhs_u, rhs_v, rhs_w
     real, dimension(:), allocatable :: xe, xc, ye, yc, ze, zc, dxe, dxc, dye, dyc, dze, dzc    ! grid edge and center coordinates for each grid-cell
     integer :: Nx, Ny, Nz, i, j, k, t_its, Vcycle
     integer, parameter :: dp = kind(1.0d0)
     real(kind=dp) :: c11, c21, c31, c12, c22, c32, c13, c23, c33
     character :: str*10
     
     ! cij, where i denote subscript and j denote superscript (i.e., RK substep no. = m)
        c11 = 0
        c21 = 0.333333333
        c31 = 0.333333333
        c12 = -0.55555556
        c22 = 0.9375
        c32 = 0.416666667
        c13 = -1.1953125
        c23 = 0.533333333
        c33 = 0.25    

     ! Lx, Ly, Lz are the lengths of the domain in the x, y, z directions
        Lx = 1
        Ly = 1
        Lz = 1

     ! Number of internal nodes    
        Nx = 128
        Ny = 128
        Nz = 128
        
     ! Reynolds Number   
        Re = 5000

     ! Current time
        time = 0

     ! Current MG cycle number
        Vcycle = 1

     ! No. of iterations in time
        t_its = 1

     ! Velocity and length scales for non-dimensionalisation
        Lo = 1
        Uo = 1   

     ! resd is the residual for L2 norm in GS iterations
        resd = 1   

     ! t_resdj is the residual for iterations in time for variable j (= u/v/w/p)
        t_resdu = 1
        t_resdv = 1
        t_resdw = 1
        t_resdp = 1

     ! t_tol is the tolerance for time iterations
        t_tol = 1.0e-6
     
     ! Kinematic viscosity
        nu = (Uo*Lo)/Re

     ! Matrix size allocation
        allocate(u(Nx+2, Ny+2, Nz+2))   ! x-velocity
        allocate(v(Nx+2, Ny+2, Nz+2))   ! y-velocity
        allocate(w(Nx+2, Ny+2, Nz+2))   ! z-velocity
        allocate(p(Nx+2, Ny+2, Nz+2))   ! pressure

        allocate(u_np1(Nx+2, Ny+2, Nz+2))   ! x-velocity at next time-step
        allocate(v_np1(Nx+2, Ny+2, Nz+2))   ! y-velocity at next time-step
        allocate(w_np1(Nx+2, Ny+2, Nz+2))   ! z-velocity at next time-step
        allocate(p_np1(Nx+2, Ny+2, Nz+2))   ! pressure at next time-step
        
        allocate(ustar(Nx+2, Ny+2, Nz+2))
        allocate(vstar(Nx+2, Ny+2, Nz+2))
        allocate(wstar(Nx+2, Ny+2, Nz+2))
        
        allocate(u1(Nx+2, Ny+2, Nz+2))
        allocate(v1(Nx+2, Ny+2, Nz+2))
        allocate(w1(Nx+2, Ny+2, Nz+2))
        allocate(p1(Nx+2, Ny+2, Nz+2))
        
        allocate(ustar2(Nx+2, Ny+2, Nz+2))
        allocate(vstar2(Nx+2, Ny+2, Nz+2))
        allocate(wstar2(Nx+2, Ny+2, Nz+2))
        
        allocate(u2(Nx+2, Ny+2, Nz+2))
        allocate(v2(Nx+2, Ny+2, Nz+2))
        allocate(w2(Nx+2, Ny+2, Nz+2))
        allocate(p2(Nx+2, Ny+2, Nz+2))
        
        allocate(ustar3(Nx+2, Ny+2, Nz+2))
        allocate(vstar3(Nx+2, Ny+2, Nz+2))
        allocate(wstar3(Nx+2, Ny+2, Nz+2))
        
        allocate(qu0(Nx+2, Ny+2, Nz+2))
        allocate(qu1(Nx+2, Ny+2, Nz+2))
        allocate(qu2(Nx+2, Ny+2, Nz+2))
        allocate(qu3(Nx+2, Ny+2, Nz+2))
        allocate(qv0(Nx+2, Ny+2, Nz+2))
        allocate(qv1(Nx+2, Ny+2, Nz+2))
        allocate(qv2(Nx+2, Ny+2, Nz+2))
        allocate(qv3(Nx+2, Ny+2, Nz+2))
        allocate(qw0(Nx+2, Ny+2, Nz+2))
        allocate(qw1(Nx+2, Ny+2, Nz+2))
        allocate(qw2(Nx+2, Ny+2, Nz+2))
        allocate(qw3(Nx+2, Ny+2, Nz+2))

        allocate(rhs_u(Nx+2, Ny+2, Nz+2))
        allocate(rhs_v(Nx+2, Ny+2, Nz+2))
        allocate(rhs_w(Nx+2, Ny+2, Nz+2))
        
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
        p = 0
        u_np1 = 0
        v_np1 = 0
        w_np1 = 0
        p_np1 = 0
        ustar = 0
        vstar = 0
        wstar = 0
        u1 = 0
        v1 = 0
        w1 = 0
        p1 = 0
        ustar2 = 0
        vstar2 = 0
        wstar2 = 0
        u2 = 0
        v2 = 0
        w2 = 0
        p2 = 0
        ustar3 = 0
        vstar3 = 0
        wstar3 = 0
        qu0 = 0
        qu1 = 0
        qu2 = 0 
        qu3 = 0
        qv0 = 0
        qv1 = 0
        qv2 = 0
        qv3 = 0
        qw0 = 0
        qw1 = 0
        qw2 = 0
        qw3 = 0
        rhs_u = 0
        rhs_v = 0
        rhs_w = 0
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

     ! Boundary Conditions
        
        ! x = 0 (Inlet)
        u_np1(1,:,:) = 0
        v_np1(1,:,:) = -v_np1(2,:,:)
        w_np1(1,:,:) = -w_np1(2,:,:)
        p_np1(1,:,:) = p_np1(2,:,:)

        ! x = Lx (Outlet)
        u_np1(Nx+1,:,:) = 0
        v_np1(Nx+2,:,:) = -v_np1(Nx+1,:,:)
        w_np1(Nx+2,:,:) = -w_np1(Nx+1,:,:)
        p_np1(Nx+2,:,:) = p_np1(Nx+1,:,:)

        ! z = 0 (ground)
        w_np1(:,:,1) = 0
        u_np1(:,:,1) = -u_np1(:,:,2)
        v_np1(:,:,1) = -v_np1(:,:,2)
        p_np1(:,:,1) = p_np1(:,:,2)

        ! z = Lz (top)
        w_np1(:,:,Nz+1) = 0
        u_np1(:,:,Nz+2) = 2 - u_np1(:,:,Nz+1)
        v_np1(:,:,Nz+2) = -v_np1(:,:,Nz+1)
        p_np1(:,:,Nz+2) = p_np1(:,:,Nz+1)

        ! y (Periodic)
        v_np1(:,1,:)    = v_np1(:,Ny,:)
        v_np1(:,Ny+1,:) = v_np1(:,2,:)
        u_np1(:,1,:)    = u_np1(:,Ny+1,:)
        u_np1(:,Ny+2,:) = u_np1(:,2,:)
        w_np1(:,1,:)    = w_np1(:,Ny+1,:)
        w_np1(:,Ny+2,:) = w_np1(:,2,:)
        p_np1(:,1,:)    = p_np1(:,Ny+1,:)
        p_np1(:,Ny+2,:) = p_np1(:,2,:)


        ! Read the input values from the last run time
!                open(unit = 600, file = "RK3CN128_Re5000_finalTime_in.dat")
!                read(600,*) time
!                do k = 1,Nz+2,1
!                        do j = 1,Ny+2,1
!                                do i = Nx+2,1
!                                           read(600,*) u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k)
!                                end do
!                        end do
!                end do
!                close(600)


!        do while ( (max(t_resdu, t_resdv, t_resdw, t_resdp) > t_tol) .and. (t_its < 1000000) )
         do while (time < 100) 
             
                ! del_t calculation
                call delta_t(del_t)

                ! update present time
                time = time + del_t

                write(*,*) 'time=',time                

                ! RK substep-1
                        call q_calc(qu1,qv1,qw1,c11,qu0,qv0,qw0,u,v,w)
         
                        ! Calculate rhs to apply thomas algorithm
                        rhs_u = u + c21*qu1
                        rhs_v = v + c21*qv1
                        rhs_w = w + c21*qw1

                        ! ustar calculation
                        call Thomas(ustar,vstar,wstar,c31,rhs_u,rhs_v,rhs_w)        

                        ! u1,v1,w1,p1 calculation
                        resd = 1
                        Vcycle = 1
                        call Multigrid(p1,u1,v1,w1,c31,ustar,vstar,wstar)

                        write(*,*) 'RK substep-1 over!'

                ! RK substep-2
                        call q_calc(qu2,qv2,qw2,c12,qu1,qv1,qw1,u1,v1,w1)

                        ! Calculate rhs to apply thomas algorithm
                        rhs_u = u1 + c22*qu2
                        rhs_v = v1 + c22*qv2
                        rhs_w = w1 + c22*qw2

                        ! ustar2 calculation
                        call Thomas(ustar2,vstar2,wstar2,c32,rhs_u,rhs_v,rhs_w)

                        ! u2,v2,w2,p2 calculation
                        resd = 1
                        Vcycle = 1
                        call Multigrid(p2,u2,v2,w2,c32,ustar2,vstar2,wstar2)

                        write(*,*) 'RK substep-2 over!'

                ! RK substep-3
                        call q_calc(qu3,qv3,qw3,c13,qu2,qv2,qw2,u2,v2,w2)

                        ! Calculate rhs to apply thomas algorithm
                        rhs_u = u2 + c23*qu3
                        rhs_v = v2 + c23*qv3
                        rhs_w = w2 + c23*qw3

                        ! ustar3 calculation
                        call Thomas(ustar3,vstar3,wstar3,c33,rhs_u,rhs_v,rhs_w)

                        ! u_np1, v_np1, w_np1, p_np1 calculation
                        resd = 1
                        Vcycle = 1
                        call Multigrid(p_np1,u_np1,v_np1,w_np1,c33,ustar3,vstar3,wstar3)

                        write(*,*) 'RK substep-3 over!'

                ! t_resd calculation
                t_resdu = 0
                t_resdv = 0
                t_resdw = 0
                t_resdp = 0

                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1
                                        t_resdu = t_resdu + (u_np1(i,j,k) - u(i,j,k))**2
                                        t_resdv = t_resdv + (v_np1(i,j,k) - v(i,j,k))**2
                                        t_resdw = t_resdw + (w_np1(i,j,k) - w(i,j,k))**2 
                                        t_resdp = t_resdp + (p_np1(i,j,k) - p(i,j,k))**2                                       
                                end do
                        end do
                end do

                t_resdu = sqrt(t_resdu/(Nx*Ny*Nz))
                t_resdv = sqrt(t_resdv/(Nx*Ny*Nz))
                t_resdw = sqrt(t_resdw/(Nx*Ny*Nz))
                t_resdp = sqrt(t_resdp/(Nx*Ny*Nz))

                write(*,*) 'time its=',t_its,'time=',time,'u res=',t_resdu,'v res=',t_resdv,'w res=',t_resdw,'p res=',t_resdp 

                if (t_its == 10000) then
                        write(*,*) 'Oops! Steady-state not reached!'
                end if      
                
                ! Updating velocity and pressure for next time-iteration
                u = u_np1
                v = v_np1
                w = w_np1
                p = p_np1
                t_its = t_its + 1 
                
        
                ! Write to output files
                ! u velocity
                open(unit = 400, file = 'RK3CN128_Re5000uvel.txt')
                do k = 1,Nz+2,1                    
                        j= (Ny+2)/2
                        write(400,*) u(:,j,k)
                end do
                close(400)

                ! w velocity
                open(unit = 300, file = 'RK3CN128_Re5000wvel.txt')
                do k = 1,Nz+2,1
                        j= (Ny+2)/2
                        write(300,*) w(:,j,k)
                end do
                close(300)
                
      
                open(unit = 500, file = "RK3CN128_Re5000_finalTime_out.dat")
                write(500,*) time
                do k = 1,Nz+2,1
                        do j = 1,Ny+2,1
                                do i = 1,Nx+2,1
                                           write(500,*) u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k)
                                end do
                        end do
                end do
                close(500)

        end do

contains

        subroutine Thomas(x1,x2,x3,c3,d1,d2,d3)
        ! x1,x2,x3 are the unknown (here, x1=u*/u**/u***, x2=v*/v**/v***, x3=w*/w**/w***)
        ! c3 is the coefficient used to calculate del_t_rkSub-step from del_t
        ! d1,d2,d3 are the rhs (explicit part) for x1(=u*), x2(=v*), x3(=w*) respectively (Ax=d)
        implicit none
                real, dimension(Nx+2,Ny+2,Nz+2), intent(in) :: d1,d2,d3
                real, dimension(Nx+2,Ny+2,Nz+2), intent(out) :: x1,x2,x3
                real :: delt_rk
                real, dimension(Nz+2) :: a, b, c, d
                integer :: i,j,k
                real(kind=dp), intent(in) :: c3

                delt_rk = c3*del_t
                
                a(1)    = 0
                c(Nz+2) = 0

                ! ustar
                do j=2,Ny+1,1                   
                        do i=2,Nx,1             
                                b(1) = 1        ! since, u(:,:,1) = -u(:,:,2)
                                c(1) = 1 
                                d(1) = 0

                                do k=2,Nz+1,1
                                        d(k) = d1(i,j,k)
                                        a(k) = (-nu*delt_rk)/(dzc(k-1)*dze(k))
                                        b(k) = 1 + (nu*delt_rk)*(1/dze(k))*( (1/dzc(k-1)) + (1/dzc(k)) )
                                        c(k) = (-nu*delt_rk)/(dzc(k)*dze(k))
                                end do

                                a(Nz+2) = 1     ! since, u(:,:,Nz+2) = 2 - u(:,:,Nz+1)
                                b(Nz+2) = 1
                                d(Nz+2) = 2
                                
                                ! Forward Elimination
                                c(1) = c(1)/b(1)
                                d(1) = d(1)/b(1)
                                do k = 2,Nz+2,1
                                        c(k) = c(k)/(b(k) - a(k)*c(k-1))
                                        d(k) = (d(k) - a(k)*d(k-1))/(b(k) - a(k)*c(k-1))
                                end do
                                
                                ! Back Substitution
                                x1(i,j,Nz+2) = d(Nz+2)
                                do k = Nz+1,1,-1
                                        x1(i,j,k) = d(k) - x1(i,j,k+1)*c(k)
                                end do

                        end do
                end do
                
                ! vstar
                do j=2,Ny,1                     
                        do i=2,Nx+1,1           
                                b(1) = 1        ! since, v(:,:,1) = -v(:,:,2)
                                c(1) = 1
                                d(1) = 0

                                do k=2,Nz+1,1
                                        d(k) = d2(i,j,k)
                                        a(k) = (-nu*delt_rk)/(dzc(k-1)*dze(k))
                                        b(k) = 1 + (nu*delt_rk)*(1/dze(k))*( (1/dzc(k-1)) + (1/dzc(k)) )
                                        c(k) = (-nu*delt_rk)/(dzc(k)*dze(k))
                                end do

                                a(Nz+2) = 1     ! since, v(:,:,Nz+2) = -v(:,:,Nz+1)
                                b(Nz+2) = 1
                                d(Nz+2) = 0

                                ! Forward Elimination
                                c(1) = c(1)/b(1)
                                d(1) = d(1)/b(1)
                                do k = 2,Nz+2,1
                                        c(k) = c(k)/(b(k) - a(k)*c(k-1))
                                        d(k) = (d(k) - a(k)*d(k-1))/(b(k) - a(k)*c(k-1))
                                end do

                                ! Back Substitution
                                x2(i,j,Nz+2) = d(Nz+2)
                                do k = Nz+1,1,-1
                                        x2(i,j,k) = d(k) - x2(i,j,k+1)*c(k)
                                end do
                        end do
                end do

                ! wstar
                do j=2,Ny+1,1                   
                        do i=2,Nx+1,1           
                                b(1) = 1        ! since, w(:,:,1) = 0
                                c(1) = 0
                                d(1) = 0

                                do k=2,Nz,1
                                        d(k) = d3(i,j,k)
                                        a(k) = (-nu*delt_rk)/(dzc(k)*dze(k))
                                        b(k) = 1 + (nu*delt_rk)*(1/dzc(k))*( (1/dze(k+1)) + (1/dze(k)) )
                                        c(k) = (-nu*delt_rk)/(dzc(k)*dze(k+1))
                                end do

                                a(Nz+1) = 0     ! since, w(:,:,Nz+1) = 0
                                b(Nz+1) = 1
                                c(Nz+1) = 0
                                d(Nz+1) = 0

                                ! Forward Elimination
                                c(1) = c(1)/b(1)
                                d(1) = d(1)/b(1)
                                do k = 2,Nz+1,1
                                        c(k) = c(k)/(b(k) - a(k)*c(k-1))
                                        d(k) = (d(k) - a(k)*d(k-1))/(b(k) - a(k)*c(k-1))
                                end do

                                ! Back Substitution
                                x3(i,j,Nz+1) = d(Nz+1)
                                do k = Nz,1,-1
                                        x3(i,j,k) = d(k) - x3(i,j,k+1)*c(k)
                                end do
                        end do
                end do

                ! Boundary conditions

                   ! x = 0 (Inlet)
                        x1(1,:,:) = 0
                        x2(1,:,:) = -x2(2,:,:)
                        x3(1,:,:) = -x3(2,:,:)

                   ! x = Lx (Outlet)
                        x1(Nx+1,:,:) = 0
                        x2(Nx+2,:,:) = -x2(Nx+1,:,:)
                        x3(Nx+2,:,:) = -x3(Nx+1,:,:)
        
                   ! z = 0 (ground)
                        x3(:,:,1) = 0
                        x1(:,:,1) = -x1(:,:,2)
                        x2(:,:,1) = -x2(:,:,2)

                   ! z = Lz (top)
                        x3(:,:,Nz+1) = 0
                        x1(:,:,Nz+2) = 2 - x1(:,:,Nz+1)
                        x2(:,:,Nz+2) = -x2(:,:,Nz+1)

                   ! y (Periodic)
                        x2(:,1,:)    = x2(:,Ny,:)
                        x2(:,Ny+1,:) = x2(:,2,:)
                        x1(:,1,:)    = x1(:,Ny+1,:)
                        x1(:,Ny+2,:) = x1(:,2,:)
                        x3(:,1,:)    = x3(:,Ny+1,:)
                        x3(:,Ny+2,:) = x3(:,2,:)

        end subroutine Thomas

        subroutine q_calc(qu_1,qv_1,qw_1,c1,qu_0,qv_0,qw_0,ut,vt,wt)

                implicit none
                real, dimension(Nx+2,Ny+2,Nz+2), intent(out) :: qu_1, qv_1, qw_1 
                real, dimension(Nx+2,Ny+2,Nz+2), intent(in) :: qu_0, qv_0, qw_0, ut, vt, wt
                real, dimension(Nx+2,Ny+2,Nz+2) :: Eu, Ev, Ew
                real :: uavg1, uavg2, vavg1, vavg2, wavg1, wavg2, d2udx2, d2udy2 
                real(kind=dp), intent(in) :: c1
                integer :: i, j, k
                
                do k = 2,Nz+1,1
                        do j = 2,Ny+1,1
                                do i = 2,Nx+1,1

                                        ! qu_1 calculation
                                        ! Convection terms
                                        uavg2 = 0.5*( ut(i+1,j,k) + ut(i,j,k) )
                                        uavg1 = 0.5*( ut(i,j,k)   + ut(i-1,j,k) )

                                        Eu(i,j,k) = -( (uavg2**2 - uavg1**2)/(dxc(i)) )

                                        vavg2 = 0.5*( vt(i+1,j,k)   + vt(i,j,k) )
                                        vavg1 = 0.5*( vt(i+1,j-1,k) + vt(i,j-1,k) )
                                        uavg2 = 0.5*( ut(i,j,k)     + ut(i,j+1,k) )
                                        uavg1 = 0.5*( ut(i,j,k)     + ut(i,j-1,k) )

                                        Eu(i,j,k) = Eu(i,j,k) - ( (uavg2*vavg2 - uavg1*vavg1)/dye(j) )

                                        wavg2 = 0.5*( wt(i+1,j,k)   + wt(i,j,k) )
                                        wavg1 = 0.5*( wt(i+1,j,k-1) + wt(i,j,k-1) )
                                        uavg2 = 0.5*( ut(i,j,k)     + ut(i,j,k+1) )
                                        uavg1 = 0.5*( ut(i,j,k)     + ut(i,j,k-1) )

                                        Eu(i,j,k) = Eu(i,j,k) - ( (uavg2*wavg2 - uavg1*wavg1)/dze(k) )
                                        
                                        ! Diffusion terms (x,y)
                                        d2udx2 = (1/dxc(i))*( (ut(i+1,j,k) - ut(i,j,k))/dxe(i+1) - &
                                                (ut(i,j,k) - ut(i-1,j,k))/dxe(i) )
                                        Eu(i,j,k) = Eu(i,j,k) + nu*d2udx2

                                        d2udy2 = (1/dye(j))*( (ut(i,j+1,k) - ut(i,j,k))/dyc(j) - &
                                                (ut(i,j,k) - ut(i,j-1,k))/dyc(j-1) )
                                        Eu(i,j,k) = Eu(i,j,k) + nu*d2udy2

                                        qu_1(i,j,k) = c1*qu_0(i,j,k) + del_t*Eu(i,j,k)

                                        ! qv_1 calculation
                                        ! Convection terms
                                        vavg2 = 0.5*( vt(i,j+1,k) + vt(i,j,k) )
                                        vavg1 = 0.5*( vt(i,j,k)   + vt(i,j-1,k) )

                                        Ev(i,j,k) = -( (vavg2**2 - vavg1**2)/(dyc(j)) )

                                        vavg2 = 0.5*( vt(i+1,j,k)   + vt(i,j,k) )
                                        vavg1 = 0.5*( vt(i-1,j,k)   + vt(i,j,k) )
                                        uavg2 = 0.5*( ut(i,j,k)     + ut(i,j+1,k) )
                                        uavg1 = 0.5*( ut(i-1,j,k)   + ut(i-1,j+1,k) )

                                        Ev(i,j,k) = Ev(i,j,k) - ( (uavg2*vavg2 - uavg1*vavg1)/dxe(i) )

                                        wavg2 = 0.5*( wt(i,j+1,k)   + wt(i,j,k) )
                                        wavg1 = 0.5*( wt(i,j+1,k-1) + wt(i,j,k-1) )
                                        vavg2 = 0.5*( vt(i,j,k)     + vt(i,j,k+1) )
                                        vavg1 = 0.5*( vt(i,j,k)     + vt(i,j,k-1) )

                                        Ev(i,j,k) = Ev(i,j,k) - ( (vavg2*wavg2 - vavg1*wavg1)/dze(k) )

                                        ! Diffusion terms (x,y)
                                        d2udx2 = (1/dxe(i))*( (vt(i+1,j,k) - vt(i,j,k))/dxc(i) - &
                                                (vt(i,j,k) - vt(i-1,j,k))/dxc(i-1) )
                                        Ev(i,j,k) = Ev(i,j,k) + nu*d2udx2

                                        d2udy2 = (1/dyc(j))*( (vt(i,j+1,k) - vt(i,j,k))/dye(j+1) - &
                                                (vt(i,j,k) - vt(i,j-1,k))/dye(j) )
                                        Ev(i,j,k) = Ev(i,j,k) + nu*d2udy2

                                        qv_1(i,j,k) = c1*qv_0(i,j,k) + del_t*Ev(i,j,k)

                                        ! qw_1 calculation
                                        ! Convection terms
                                        wavg2 = 0.5*( wt(i,j,k+1) + wt(i,j,k) )
                                        wavg1 = 0.5*( wt(i,j,k)   + wt(i,j,k-1) )

                                        Ew(i,j,k) = -( (wavg2**2 - wavg1**2)/(dzc(k)) )

                                        wavg2 = 0.5*( wt(i+1,j,k)   + wt(i,j,k) )
                                        wavg1 = 0.5*( wt(i-1,j,k)   + wt(i,j,k) )
                                        uavg2 = 0.5*( ut(i,j,k)     + ut(i,j,k+1) )
                                        uavg1 = 0.5*( ut(i-1,j,k)   + ut(i-1,j,k+1) )

                                        Ew(i,j,k) = Ew(i,j,k) - ( (uavg2*wavg2 - uavg1*wavg1)/dxe(i) )

                                        wavg2 = 0.5*( wt(i,j+1,k)   + wt(i,j,k) )
                                        wavg1 = 0.5*( wt(i,j-1,k)   + wt(i,j,k) )
                                        vavg2 = 0.5*( vt(i,j,k)     + vt(i,j,k+1) )
                                        vavg1 = 0.5*( vt(i,j-1,k)   + vt(i,j-1,k+1) )

                                        Ew(i,j,k) = Ew(i,j,k) - ( (vavg2*wavg2 - vavg1*wavg1)/dye(j) )

                                        ! Diffusion terms (x,y)
                                        d2udx2 = (1/dxe(i))*( (wt(i+1,j,k) - wt(i,j,k))/dxc(i) - &
                                                 (wt(i,j,k) - wt(i-1,j,k))/dxc(i-1) )
                                        Ew(i,j,k) = Ew(i,j,k) + nu*d2udx2

                                        d2udy2 = (1/dye(j))*( (wt(i,j+1,k) - wt(i,j,k))/dyc(j) - &
                                                 (wt(i,j,k) - wt(i,j-1,k))/dyc(j-1) )
                                        Ew(i,j,k) = Ew(i,j,k) + nu*d2udy2

                                        qw_1(i,j,k) = c1*qw_0(i,j,k) + del_t*Ew(i,j,k)

                                end do
                        end do
                end do

        end subroutine q_calc

        subroutine Multigrid(pt,ut,vt,wt,c3,us,vs,ws)
                ! pt, ut, vt, wt can be p1/p2/p, u1/u2/u, v1/v2/v, w1/w2/w respectively
                ! c3 is the coefficient used to calculate del_t_rkSub-step from del_t
                ! us, vs, ws can be ustar/ustar2/ustar3, vstar/vstar2/vstar3, wstar/wstar2/wstar3 respectively
                
                implicit none
                real, dimension(Nx+2,Ny+2,Nz+2), intent(out) :: pt, ut, vt, wt
                real, dimension(Nx+2,Ny+2,Nz+2), intent(in) :: us, vs, ws
                real(kind=dp), intent(in) :: c3
                real :: delt_rk

                ! Going upto 5 grid levels
                real, dimension(:,:,:), allocatable :: phi1, rhs1, res1, phi2, rhs2, res2, phi3,rhs3, res3, phi4, rhs4, res4, &
                                                        phi5, rhs5, res5, phi54, phi43, phi32, phi21
                real, dimension(:), allocatable :: dxe2, dxc2, dye2, dyc2, dze2, dzc2, dxe3, dxc3, dye3, dyc3, dze3, dzc3, &
                                                        dxe4, dxc4, dye4, dyc4, dze4, dzc4, dxe5, dxc5, dye5, dyc5, dze5, dzc5
                ! dxcj, dxej (where j = 2/3/4/5) are the distances b/w cell centers and cell edges for different MG levels 
                real :: tol
                integer :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4, nx5, ny5, nz5
                integer :: i,j,k, plot_var, ic, jc, kc, i_f, jf, kf

                delt_rk = c3*del_t
                
                ! Tolerance
                tol = 5.0e-4

                ! No. of V cycles traversed
                Vcycle = 1

                ! nx1, ny1, nz1 are the number of internal nodes in the finest mesh (level-1)
                ! nx2, ny2, nz2 are the number of internal nodes at level-2
                ! nx3, ny3, nz3 are the number of internal nodes at level-3
                ! nx4, ny4, nz4 are the number of internal nodes at level-4
                ! nx5, ny5, nz5 are the number of internal nodes at level-5
                nx1 = Nx
                ny1 = Ny
                nz1 = Nz

                nx2 = nx1/2
                ny2 = ny1/2
                nz2 = nz1/2

                nx3 = nx2/2
                ny3 = ny2/2
                nz3 = nz2/2

                nx4 = nx3/2
                ny4 = ny3/2
                nz4 = nz3/2

                nx5 = nx4/2
                ny5 = ny4/2
                nz5 = nz4/2

                ! Matrix size allocation
                allocate(phi1(nx1+2, ny1+2, nz1+2))
                allocate(rhs1(nx1+2, ny1+2, nz1+2))
                allocate(res1(nx1+2, ny1+2, nz1+2))
                allocate(phi21(nx1+2, ny1+2, nz1+2))

                allocate(phi2(nx2+2, ny2+2, nz2+2))
                allocate(rhs2(nx2+2, ny2+2, nz2+2))
                allocate(res2(nx2+2, ny2+2, nz2+2))
                allocate(phi32(nx2+2, ny2+2, nz2+2))

                allocate(phi3(nx3+2, ny3+2, nz3+2))
                allocate(rhs3(nx3+2, ny3+2, nz3+2))
                allocate(res3(nx3+2, ny3+2, nz3+2))
                allocate(phi43(nx3+2, ny3+2, nz3+2))

                allocate(phi4(nx4+2, ny4+2, nz4+2))
                allocate(rhs4(nx4+2, ny4+2, nz4+2))
                allocate(res4(nx4+2, ny4+2, nz4+2))
                allocate(phi54(nx4+2, ny4+2, nz4+2))

                allocate(phi5(nx5+2, ny5+2, nz5+2))
                allocate(rhs5(nx5+2, ny5+2, nz5+2))
                allocate(res5(nx5+2, ny5+2, nz5+2))

                allocate(dxe2(nx2+2))
                allocate(dxc2(nx2+2))
                allocate(dye2(ny2+2))
                allocate(dyc2(ny2+2))
                allocate(dze2(nz2+2))
                allocate(dzc2(nz2+2))

                allocate(dxe3(nx3+2))
                allocate(dxc3(nx3+2))
                allocate(dye3(ny3+2))
                allocate(dyc3(ny3+2))
                allocate(dze3(nz3+2))
                allocate(dzc3(nz3+2))
                
                allocate(dxe4(nx4+2))
                allocate(dxc4(nx4+2))
                allocate(dye4(ny4+2))
                allocate(dyc4(ny4+2))
                allocate(dze4(nz4+2))
                allocate(dzc4(nz4+2))

                allocate(dxe5(nx5+2))
                allocate(dxc5(nx5+2))
                allocate(dye5(ny5+2))
                allocate(dyc5(ny5+2))
                allocate(dze5(nz5+2))
                allocate(dzc5(nz5+2))

                ! Initialisation
                phi1 = 0
                rhs1 = 0
                res1 = 0
                phi21 = 0

                phi2 = 0
                rhs2 = 0
                res2 = 0
                phi32 = 0

                phi3 = 0
                rhs3 = 0
                res3 = 0
                phi43 = 0

                phi4 = 0
                rhs4 = 0
                res4 = 0
                phi54 = 0

                phi5 = 0
                rhs5 = 0
                res5 = 0

                dxe2 = 0
                dxc2 = 0
                dye2 = 0
                dyc2 = 0
                dze2 = 0
                dzc2 = 0

                dxe3 = 0
                dxc3 = 0
                dye3 = 0
                dyc3 = 0
                dze3 = 0
                dzc3 = 0

                dxe4 = 0
                dxc4 = 0
                dye4 = 0
                dyc4 = 0
                dze4 = 0
                dzc4 = 0

                dxe5 = 0
                dxc5 = 0
                dye5 = 0
                dyc5 = 0
                dze5 = 0
                dzc5 = 0
                
                ! Grid-size calculation for lower levels (2 to 5)
                
                ! x-direction
                ! level-2
                do ic = 1,nx2+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dxc2(ic) = dxe(2) + dxc(2) 
                        else if (ic == nx2+1) then
                              dxc2(ic) = dxc(nx1) + dxe(nx1+1)
                        else
                              dxc2(ic) = dxc(i_f)/2 + dxc(i_f+1) + dxc(i_f+2)/2
                        end if

                        dxe2(ic+1) = dxc2(ic)
                end do

                ! level-3
                do ic = 1,nx3+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dxc3(ic) = dxe2(2) + dxc2(2)
                        else if (ic == nx3+1) then
                              dxc3(ic) = dxc2(nx2) + dxe2(nx2+1)
                        else
                              dxc3(ic) = dxc2(i_f)/2 + dxc2(i_f+1) + dxc2(i_f+2)/2
                        end if

                        dxe3(ic+1) = dxc3(ic)
                end do

                ! level-4
                do ic = 1,nx4+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dxc4(ic) = dxe3(2) + dxc3(2)
                        else if (ic == nx4+1) then
                              dxc4(ic) = dxc3(nx3) + dxe3(nx3+1)
                        else
                              dxc4(ic) = dxc3(i_f)/2 + dxc3(i_f+1) + dxc3(i_f+2)/2
                        end if

                        dxe4(ic+1) = dxc4(ic)
                end do

                ! level-5
                do ic = 1,nx5+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dxc5(ic) = dxe4(2) + dxc4(2)
                        else if (ic == nx5+1) then
                              dxc5(ic) = dxc4(nx4) + dxe4(nx4+1)
                        else
                              dxc5(ic) = dxc4(i_f)/2 + dxc4(i_f+1) + dxc4(i_f+2)/2
                        end if

                        dxe5(ic+1) = dxc5(ic)
                end do

                ! y-direction
                ! level-2
                do ic = 1,ny2+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dyc2(ic) = dye(2) + dyc(2)
                        else if (ic == ny2+1) then
                              dyc2(ic) = dyc(ny1) + dye(ny1+1)
                        else
                              dyc2(ic) = dyc(i_f)/2 + dyc(i_f+1) + dyc(i_f+2)/2
                        end if

                        dye2(ic+1) = dyc2(ic)
                end do

                ! level-3
                do ic = 1,ny3+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dyc3(ic) = dye2(2) + dyc2(2)
                        else if (ic == ny3+1) then
                              dyc3(ic) = dyc2(ny2) + dye2(ny2+1)
                        else
                              dyc3(ic) = dyc2(i_f)/2 + dyc2(i_f+1) + dyc2(i_f+2)/2
                        end if

                        dye3(ic+1) = dyc3(ic)
                end do

                ! level-4
                do ic = 1,ny4+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dyc4(ic) = dye3(2) + dyc3(2)
                        else if (ic == ny4+1) then
                              dyc4(ic) = dyc3(ny3) + dye3(ny3+1)
                        else
                              dyc4(ic) = dyc3(i_f)/2 + dyc3(i_f+1) + dyc3(i_f+2)/2
                        end if

                        dye4(ic+1) = dyc4(ic)
                end do

                ! level-5
                do ic = 1,ny5+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dyc5(ic) = dye4(2) + dyc4(2)
                        else if (ic == ny5+1) then
                              dyc5(ic) = dyc4(ny4) + dye4(ny4+1)
                        else
                              dyc5(ic) = dyc4(i_f)/2 + dyc4(i_f+1) + dyc4(i_f+2)/2
                        end if

                        dye5(ic+1) = dyc5(ic)
                end do

                ! z-direction
                ! level-2
                do ic = 1,nz2+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dzc2(ic) = dze(2) + dzc(2)
                        else if (ic == nz2+1) then
                              dzc2(ic) = dzc(nz1) + dze(nz1+1)
                        else
                              dzc2(ic) = dzc(i_f)/2 + dzc(i_f+1) + dzc(i_f+2)/2
                        end if

                        dze2(ic+1) = dzc2(ic)
                end do

                ! level-3
                do ic = 1,nz3+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dzc3(ic) = dze2(2) + dzc2(2)
                        else if (ic == nz3+1) then
                              dzc3(ic) = dzc2(nz2) + dze2(nz2+1)
                        else
                              dzc3(ic) = dzc2(i_f)/2 + dzc2(i_f+1) + dzc2(i_f+2)/2
                        end if

                        dze3(ic+1) = dzc3(ic)
                end do

                ! level-4
                do ic = 1,nz4+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dzc4(ic) = dze3(2) + dzc3(2)
                        else if (ic == nz4+1) then
                              dzc4(ic) = dzc3(nz3) + dze3(nz3+1)
                        else
                              dzc4(ic) = dzc3(i_f)/2 + dzc3(i_f+1) + dzc3(i_f+2)/2
                        end if

                        dze4(ic+1) = dzc4(ic)
                end do

                ! level-5
                do ic = 1,nz5+1,1
                        i_f = (ic-2) + ic
                        if (ic == 1) then
                              dzc5(ic) = dze4(2) + dzc4(2)
                        else if (ic == nz5+1) then
                              dzc5(ic) = dzc4(nz4) + dze4(nz4+1)
                        else
                              dzc5(ic) = dzc4(i_f)/2 + dzc4(i_f+1) + dzc4(i_f+2)/2
                        end if

                        dze5(ic+1) = dzc5(ic)
                end do


                ! rhs1 calculation
                do k = 2,nz1+1,1
                        do j = 2,ny1+1,1
                                do i = 2,nx1+1,1
                                        rhs1(i,j,k) = (1/delt_rk)*( ((us(i,j,k)-us(i-1,j,k))/dxe(i)) + &
                                                      ((vs(i,j,k)-vs(i,j-1,k))/dye(j)) + ((ws(i,j,k)-ws(i,j,k-1))/dze(k)) )
                                end do
                        end do
                end do

             
                resd = 1

                Do while (resd > tol)
                        write(200,*) 'V cycle=',Vcycle,':'
                        write(*,*) 'V cycle=', Vcycle, ':'        

                        call GS(1,nx1,ny1,nz1,6,dxc,dxe,dyc,dye,dzc,dze,phi1,rhs1,res1)
                        
                        call restrict(res1,nx1,ny1,nz1,rhs2,nx2,ny2,nz2)
                        call GS(2,nx2,ny2,nz2,5,dxc2,dxe2,dyc2,dye2,dzc2,dze2,phi2,rhs2,res2)
                        
                        call restrict(res2,nx2,ny2,nz2,rhs3,nx3,ny3,nz3)
                        call GS(3,nx3,ny3,nz3,5,dxc3,dxe3,dyc3,dye3,dzc3,dze3,phi3,rhs3,res3)
                        
                        call restrict(res3,nx3,ny3,nz3,rhs4,nx4,ny4,nz4)
                        call GS(4,nx4,ny4,nz4,4,dxc4,dxe4,dyc4,dye4,dzc4,dze4,phi4,rhs4,res4)
                        
                        call restrict(res4,nx4,ny4,nz4,rhs5,nx5,ny5,nz5)
                        call GS(5,nx5,ny5,nz5,4,dxc5,dxe5,dyc5,dye5,dzc5,dze5,phi5,rhs5,res5)
                        
                        call prolongate(phi54,nx4,ny4,nz4,phi5,nx5,ny5,nz5)
                        phi4 = phi4 + phi54
                        call GS(4,nx4,ny4,nz4,4,dxc4,dxe4,dyc4,dye4,dzc4,dze4,phi4,rhs4,res4)

                        call prolongate(phi43,nx3,ny3,nz3,phi4,nx4,ny4,nz4)
                        phi3 = phi3 + phi43
                        call GS(3,nx3,ny3,nz3,5,dxc3,dxe3,dyc3,dye3,dzc3,dze3,phi3,rhs3,res3)

                        call prolongate(phi32,nx2,ny2,nz2,phi3,nx3,ny3,nz3)
                        phi2 = phi2 + phi32
                        call GS(2,nx2,ny2,nz2,5,dxc2,dxe2,dyc2,dye2,dzc2,dze2,phi2,rhs2,res2)

                        call prolongate(phi21,nx1,ny1,nz1,phi2,nx2,ny2,nz2)
                        phi1 = phi1 + phi21
                        call GS(1,nx1,ny1,nz1,6,dxc,dxe,dyc,dye,dzc,dze,phi1,rhs1,res1)

                        Vcycle = Vcycle + 1
                 end do
                
                pt = phi1
                
                ! ut calculation
                do k = 2,nz1+1,1
                        do j = 2,ny1+1,1
                                do i = 2,nx1,1
                                        ut(i,j,k) = us(i,j,k) - (delt_rk/dxc(i))*(pt(i+1,j,k) - pt(i,j,k))          
                                end do
                        end do
                end do

                ! vt calculation
                do k = 2,nz1+1,1
                        do j = 2,ny1,1
                                do i = 2,nx1+1,1
                                        vt(i,j,k) = vs(i,j,k) - (delt_rk/dyc(j))*(pt(i,j+1,k) - pt(i,j,k))
                                end do
                        end do
                end do

                ! wt calculation
                do k = 2,nz1,1
                        do j = 2,ny1+1,1
                                do i = 2,nx1+1,1
                                        wt(i,j,k) = ws(i,j,k) - (delt_rk/dzc(k))*(pt(i,j,k+1) - pt(i,j,k))
                                end do
                        end do
                end do

                ! Boundary Conditions for ut,vt,wt
                ! x = 0 (Inlet)
                ut(1,:,:) = 0
                vt(1,:,:) = -vt(2,:,:)
                wt(1,:,:) = -wt(2,:,:)

                ! x = Lx (Outlet)
                ut(nx1+1,:,:) = 0
                vt(nx1+2,:,:) = -vt(nx1+1,:,:)
                wt(nx1+2,:,:) = -wt(nx1+1,:,:)

                ! z = 0 (ground)
                wt(:,:,1) = 0
                ut(:,:,1) = -ut(:,:,2)
                vt(:,:,1) = -vt(:,:,2)

                ! z = Lz (top)
                wt(:,:,nz1+1) = 0
                ut(:,:,nz1+2) = 2 - ut(:,:,nz1+1)
                vt(:,:,nz1+2) = -vt(:,:,nz1+1)

                ! y (Periodic)
                vt(:,1,:)     = vt(:,ny1,:)
                vt(:,ny1+1,:) = vt(:,2,:)
                ut(:,1,:)     = ut(:,ny1+1,:)
                ut(:,ny1+2,:) = ut(:,2,:)
                wt(:,1,:)     = wt(:,ny1+1,:)
                wt(:,ny1+2,:) = wt(:,2,:)

                deallocate(phi1)
                deallocate(rhs1)
                deallocate(res1)
                deallocate(phi21)

                deallocate(phi2)
                deallocate(rhs2)
                deallocate(res2)
                deallocate(phi32)

                deallocate(phi3)
                deallocate(rhs3)
                deallocate(res3)
                deallocate(phi43)

                deallocate(phi4)
                deallocate(rhs4)
                deallocate(res4)
                deallocate(phi54)

                deallocate(phi5)
                deallocate(rhs5)
                deallocate(res5) 

                deallocate(dxe2)
                deallocate(dxc2)
                deallocate(dye2)
                deallocate(dyc2)
                deallocate(dze2)
                deallocate(dzc2)

                deallocate(dxe3)
                deallocate(dxc3)
                deallocate(dye3)
                deallocate(dyc3)
                deallocate(dze3)
                deallocate(dzc3)

                deallocate(dxe4)
                deallocate(dxc4)
                deallocate(dye4)
                deallocate(dyc4)
                deallocate(dze4)
                deallocate(dzc4)

                deallocate(dxe5)
                deallocate(dxc5)
                deallocate(dye5)
                deallocate(dyc5)
                deallocate(dze5)
                deallocate(dzc5)

        end subroutine Multigrid

        subroutine restrict(res_l1,Nx1,Ny1,Nz1,rhs_l2,Nx2,Ny2,Nz2)

                ! Nx1, Ny1, Nz1 are the number of internal nodes of res_l1 (fine)
                ! Nx2, Ny2, Nz2 are the number of internal nodes of rhs_l2 (coarse)
                
                implicit none
                integer :: ic, jc, kc, i_f, jf, kf, Nx1, Ny1, Nz1, Nx2, Ny2, Nz2
                real, dimension(Nx2+2,Ny2+2,Nz2+2), intent(out) :: rhs_l2
                real, dimension(Nx1+2,Ny1+2,Nz1+2) :: res_l1

                Do kc = 2,Nz2+1,1
                        Do jc = 2,Ny2+1,1
                                Do ic = 2,Nx2+1,1
                                        i_f = (ic-2) + ic
                                        jf = (jc-2) + jc
                                        kf = (kc-2) + kc
                                        rhs_l2(ic,jc,kc) = 0.128*( res_l1(i_f,jf,kf) + res_l1(i_f+1,jf,kf) + &
                                                                res_l1(i_f+1,jf+1,kf) + res_l1(i_f,jf+1,kf) + &
                                                                res_l1(i_f,jf,kf+1) + res_l1(i_f+1,jf,kf+1) + &
                                                                res_l1(i_f+1,jf+1,kf+1) + res_l1(i_f,jf+1,kf+1) )
                                end do
                        end do
                end do
                
        end subroutine restrict

        subroutine prolongate(phi_21,Nx1,Ny1,Nz1,phi_2,Nx2,Ny2,Nz2)

                ! Nx1, Ny1, Nz1 are the number of internal nodes of phi_21 (fine)
                ! Nx2, Ny2, Nz2 are the number of internal nodes of phi_2 (coarse)
                
                implicit none
                integer :: ic, jc, kc, i_f, jf, kf, Nx1, Ny1, Nz1, Nx2, Ny2, Nz2
                real, dimension(Nx1+2,Ny1+2,Nz1+2), intent(out) :: phi_21
                real, dimension(Nx2+2,Ny2+2,Nz2+2) :: phi_2
                
                Do kf = 2,Nz1+2,2
                        Do jf = 2,Ny1+2,2
                                Do i_f = 2,Nx1+2,2
                                        ic = i_f/2
                                        jc = jf/2
                                        kc = kf/2
                                        phi_21(i_f-1,jf-1,kf-1) = ( 27*phi_2(ic,jc,kc) + 9*phi_2(ic+1,jc,kc) + &
                                                                9*phi_2(ic,jc,kc+1) + 9*phi_2(ic,jc+1,kc) + &
                                                                3*phi_2(ic+1,jc+1,kc) + 3*phi_2(ic+1,jc,kc+1) + &
                                                                3*phi_2(ic,jc+1,kc+1) + phi_2(ic+1,jc+1,kc+1) )/64

                                        phi_21(i_f,jf-1,kf-1) = ( 27*phi_2(ic+1,jc,kc) + 9*phi_2(ic+1,jc+1,kc) + &
                                                                9*phi_2(ic,jc,kc) + 9*phi_2(ic+1,jc,kc+1) + &
                                                                3*phi_2(ic+1,jc+1,kc+1) + 3*phi_2(ic,jc,kc+1) + &
                                                                3*phi_2(ic,jc+1,kc) + phi_2(ic,jc+1,kc+1) )/64

                                        phi_21(i_f,jf,kf-1) = ( 27*phi_2(ic+1,jc+1,kc) + 9*phi_2(ic+1,jc,kc) + &
                                                                9*phi_2(ic+1,jc+1,kc+1) + 9*phi_2(ic,jc+1,kc) + &
                                                                3*phi_2(ic+1,jc,kc+1) + 3*phi_2(ic,jc+1,kc+1) + &
                                                                3*phi_2(ic,jc,kc) + phi_2(ic,jc,kc+1) )/64

                                        phi_21(i_f-1,jf,kf-1) = ( 27*phi_2(ic,jc+1,kc) + 9*phi_2(ic,jc+1,kc+1) + &
                                                                9*phi_2(ic,jc,kc) + 9*phi_2(ic+1,jc+1,kc) + &
                                                                3*phi_2(ic+1,jc,kc) + 3*phi_2(ic+1,jc+1,kc+1) + &
                                                                3*phi_2(ic,jc,kc+1) + phi_2(ic+1,jc,kc+1) )/64

                                        phi_21(i_f-1,jf-1,kf) = ( 27*phi_2(ic,jc,kc+1) + 9*phi_2(ic+1,jc,kc+1) + &
                                                                9*phi_2(ic,jc+1,kc+1) + 9*phi_2(ic,jc,kc) + &
                                                                3*phi_2(ic,jc+1,kc) + 3*phi_2(ic+1,jc,kc) + &
                                                                3*phi_2(ic+1,jc+1,kc+1) + phi_2(ic+1,jc+1,kc) )/64

                                        phi_21(i_f,jf-1,kf) = ( 27*phi_2(ic+1,jc,kc+1) + 9*phi_2(ic+1,jc,kc) + &
                                                                9*phi_2(ic,jc,kc+1) + 9*phi_2(ic+1,jc+1,kc+1) + &
                                                                3*phi_2(ic+1,jc+1,kc) + 3*phi_2(ic,jc,kc) + &
                                                                3*phi_2(ic,jc+1,kc+1) + phi_2(ic,jc+1,kc) )/64

                                        phi_21(i_f,jf,kf) = ( 27*phi_2(ic+1,jc+1,kc+1) + 9*phi_2(ic+1,jc+1,kc) + &
                                                                9*phi_2(ic,jc+1,kc+1) + 9*phi_2(ic+1,jc,kc+1) + &
                                                                3*phi_2(ic+1,jc,kc) + 3*phi_2(ic,jc,kc+1) + &
                                                                3*phi_2(ic,jc+1,kc) + phi_2(ic,jc,kc) )/64

                                        phi_21(i_f-1,jf,kf) = ( 27*phi_2(ic,jc+1,kc+1) + 9*phi_2(ic,jc+1,kc) + &
                                                                9*phi_2(ic,jc,kc+1) + 9*phi_2(ic+1,jc+1,kc+1) + &
                                                                3*phi_2(ic+1,jc+1,kc) + 3*phi_2(ic+1,jc,kc+1) + &
                                                                3*phi_2(ic,jc,kc) + phi_2(ic+1,jc,kc) )/64
                
                                end do
                        end do
                end do
                
        end subroutine prolongate

        subroutine GS(MGlevel,nxt,nyt,nzt,nits,del_xc,del_xe,del_yc,del_ye,del_zc,del_ze,phi_old,rhs,res)
                ! nits = No. of desired GS iterations at any MG level 
                implicit none
                integer :: nxt, nyt, nzt, i, j, k, its, nits, MGlevel
                real, dimension(nxt+2,nyt+2,nzt+2) :: phi, phi_diff, c, cw, ce, cs, cn, ca, cb, rhs
                real, dimension(nxt+2,nyt+2,nzt+2), intent(out) :: phi_old, res
                real, dimension(nxt+2) :: del_xe, del_xc
                real, dimension(nyt+2) :: del_ye, del_yc
                real, dimension(nzt+2) :: del_ze, del_zc

                ! Number of internal nodes at any given MG level = nxt, nyt, nzt. Total nodes = nxt+2, nyt+2, nzt+2 

                ! Initialising all arrays to zero
                phi = 0
                phi_diff = 0
                c = 0
                cw = 0
                ce = 0
                cn = 0
                cs = 0
                ca = 0
                cb = 0

                ! Saving the coefficients
                do i = 2,nxt+1,1
                        do j = 2,nyt+1,1
                                do k = 2,nzt+1,1
                                        c(i,j,k) = -( (1/(del_xc(i)*del_xe(i)) + 1/(del_xc(i-1)*del_xe(i))) + &
                                                        (1/(del_yc(j)*del_ye(j)) + 1/(del_yc(j-1)*del_ye(j))) + &
                                                        (1/(del_zc(k)*del_ze(k)) + 1/(del_zc(k-1)*del_ze(k))) )
                                        cw(i,j,k) = 1/(del_xc(i-1)*del_xe(i))
                                        ce(i,j,k) = 1/(del_xc(i)*del_xe(i))
                                        cs(i,j,k) = 1/(del_yc(j-1)*del_ye(j))
                                        cn(i,j,k) = 1/(del_yc(j)*del_ye(j))
                                        ca(i,j,k) = 1/(del_zc(k)*del_ze(k))
                                        cb(i,j,k) = 1/(del_zc(k-1)*del_ze(k))
                                end do
                        end do
                end do

                resd = 1
                its = 0

                ! Boundary Conditions
                if (MGlevel == 1) then
                                ! x = 0 (Inlet)
                                phi(1,:,:) = phi(2,:,:)

                                ! x = Lx (Outlet)
                                phi(nxt+2,:,:) = phi(nxt+1,:,:)

                                ! z = 0 (ground)
                                phi(:,:,1) = phi(:,:,2)

                                ! z = Lz (top)
                                phi(:,:,nzt+2) = phi(:,:,nzt+1)

                                ! y (Periodic)
                                phi(:,1,:)    = phi(:,nyt+1,:)
                                phi(:,nyt+2,:) = phi(:,2,:)
                else
                                phi(1,:,:)     = 0
                                phi(nxt+2,:,:) = 0 
                                phi(:,:,1)     = 0
                                phi(:,:,nzt+2) = 0
                                phi(:,1,:)     = 0
                                phi(:,nyt+2,:) = 0
                end if


                ! GS iterations
                do while( its < nits )
                        if (MGlevel == 1) then
                                do k = 2,nzt+1,1
                                        do j = 2,nyt+1,1
                                                do i = 2,nxt+1,1
                                                        phi(i,j,k) = (rhs(i,j,k) - cw(i,j,k)*phi(i-1,j,k) - &
                                                                        ce(i,j,k)*phi_old(i+1,j,k) - &
                                                                        cs(i,j,k)*phi(i,j-1,k) - cn(i,j,k)*phi_old(i,j+1,k) - &
                                                                        cb(i,j,k)*phi(i,j,k-1) - &
                                                                        ca(i,j,k)*phi_old(i,j,k+1) )/c(i,j,k)                        
                                                end do
                                        end do     
                                end do

                                ! Boundary Conditions
                                ! x = 0 (Inlet)
                                phi(1,:,:) = phi(2,:,:)

                                ! x = Lx (Outlet)
                                phi(nxt+2,:,:) = phi(nxt+1,:,:)

                                ! z = 0 (ground)
                                phi(:,:,1) = phi(:,:,2)

                                ! z = Lz (top)
                                phi(:,:,nzt+2) = phi(:,:,nzt+1)

                                ! y (Periodic)
                                phi(:,1,:)    = phi(:,nyt+1,:)
                                phi(:,nyt+2,:) = phi(:,2,:)


                        else
                                do k = 2,nzt+1,1
                                        do j = 2,nyt+1,1
                                                do i = 2,nxt+1,1
                                                        phi(i,j,k) = (rhs(i,j,k) - cw(i,j,k)*phi(i-1,j,k) - &
                                                                        ce(i,j,k)*phi_old(i+1,j,k) - &
                                                                        cs(i,j,k)*phi(i,j-1,k) - cn(i,j,k)*phi_old(i,j+1,k) - &
                                                                        cb(i,j,k)*phi(i,j,k-1) - &
                                                                        ca(i,j,k)*phi_old(i,j,k+1) )/c(i,j,k)
                                                end do
                                        end do        
                                end do

                                ! Boundary Conditions
                                phi(1,:,:)     = 0
                                phi(nxt+2,:,:) = 0
                                phi(:,:,1)     = 0
                                phi(:,:,nzt+2) = 0
                                phi(:,1,:)     = 0
                                phi(:,nyt+2,:) = 0
                        end if                  
                
                        phi_diff = phi - phi_old
                        resd = 0

                        ! Calculating residual
                        do k = 2,nzt+1,1
                                do j = 2,nyt+1,1
                                        do i = 2,nxt+1,1
                                                resd = resd + (phi_diff(i,j,k))**2
                                        end do
                                end do
                        end do

                        resd = sqrt(resd/(nxt*nyt*nzt))
                        phi_old = phi
                        its = its + 1

                write(*,*) 'time its=',t_its,'V cycle=',Vcycle, 'MG level=',MGlevel, 'GS its=',its, 'res=',resd

                end do

                do k = 2,nzt+1,1
                        do j = 2,nyt+1,1
                                do i = 2,nxt+1,1
                                        res(i,j,k) = rhs(i,j,k) - cw(i,j,k)*phi(i-1,j,k) - ce(i,j,k)*phi(i+1,j,k) - &
                                                     cs(i,j,k)*phi(i,j-1,k) - cn(i,j,k)*phi(i,j+1,k) - cb(i,j,k)*phi(i,j,k-1) - &
                                                     ca(i,j,k)*phi(i,j,k+1) - c(i,j,k)*phi(i,j,k)
                                end do
                        end do
                end do

        end subroutine GS

        subroutine delta_t(delt)
        
                implicit none
                real :: delt, delt_x, delt_y, delt_z, delx, dely, delz
                delx = minval(dxe(2:Nx+2))
                dely = minval(dye(2:Ny+2))
                delz = minval(dze(2:Nz+2))

                delt_x = (1.2*delx)/maxval(u_np1)
                delt_y = (1.2*dely)/maxval(v_np1)
                delt_z = (1.2*delz)/maxval(w_np1)

                delt = min(delt_x, delt_y, delt_z)

        end subroutine delta_t

end program main 
