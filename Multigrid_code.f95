! Multigrid code by Shreya Agrawal (160662)
! in Fortran 95
! Problem is stated in the Mutigrid-report pdf file.

program main

implicit none

! Going upto 5 grid levels
real, dimension(:,:,:), allocatable :: phi1, rhs1, res1, phi2, rhs2, res2,&
phi3,rhs3, res3, phi4, rhs4, res4, phi5, rhs5, res5, phi54, phi43, phi32, phi21
real :: Lx, Ly, Lz, tol, del_xt, del_yt, del_zt, resd
! del_xt, del_yt, del_zt are temporary variables
integer :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4, nx5, ny5, nz5
integer :: i,j,k, Vcycle, plot_var

! Tolerance
tol = 1.0e-5

! resd is the residual for L2 norm in GS iterations
resd = 1

! Plot variable = useful in naming output files
plot_var = 1

! Lx, Ly, Lz are the lengths of the domain in the x, y, z directions
Lx = 1
Ly = 1
Lz = 1

! No. of V cycles traversed
Vcycle = 1 

! nx1, ny1, nz1 are the number of internal nodes in the finest mesh (level-1)
! nx2, ny2, nz2 are the number of internal nodes at level-2
! nx3, ny3, nz3 are the number of internal nodes at level-3
! nx4, ny4, nz4 are the number of internal nodes at level-4
! nx5, ny5, nz5 are the number of internal nodes at level-5
nx1 = 64
ny1 = 64
nz1 = 64

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

! rhs1 calculation
do k = 2,nz1+1,1
        do j = 2,ny1+1,1
                do i = 2,nx1+1,1
                        del_xt = Lx/nx1
                        del_zt = Lz/nz1
                        rhs1(i,j,k) = 50000*( (100*( (1-(del_xt*(0.5+i-2)))**2 + &
                                   (del_zt*(0.5+k-2))**2 ))-2 )*exp(-50*((1-(del_xt*(0.5+i-2)))**2 + &
                                   (del_zt*(0.5+k-2))**2 ))
                end do
        end do
end do

! fort.6 will have residual information for GS iteration
open(unit = 200, file = 'fort.6')

! V cycle iterations
Do while (resd > tol)
	call GS(1,nx1,ny1,nz1,10,phi1,rhs1,res1)

	if (Vcycle == 2) then
		call plot('MGV02L1_D.txt',phi1,nx1,ny1,nz1)

	else if (Vcycle == 30) then
		call plot('MGV30L1_D.txt',phi1,nx1,ny1,nz1)
	else
	end if
 
!	if( resd < tol ) exit
	call restrict(res1,nx1,ny1,nz1,rhs2,nx2,ny2,nz2)

	call GS(2,nx2,ny2,nz2,8,phi2,rhs2,res2)

        if (Vcycle == 2) then
                call plot('MGV02L2_D.txt',phi2,nx2,ny2,nz2)

        else if (Vcycle == 30) then
                call plot('MGV30L2_D.txt',phi2,nx2,ny2,nz2)
        else
        end if

	call restrict(res2,nx2,ny2,nz2,rhs3,nx3,ny3,nz3)

	call GS(3,nx3,ny3,nz3,6,phi3,rhs3,res3)

        if (Vcycle == 2) then
                call plot('MGV02L3_D.txt',phi3,nx3,ny3,nz3)

        else if (Vcycle == 30) then
                call plot('MGV30L3_D.txt',phi3,nx3,ny3,nz3)
        else
        end if

	call restrict(res3,nx3,ny3,nz3,rhs4,nx4,ny4,nz4)

	call GS(4,nx4,ny4,nz4,4,phi4,rhs4,res4)

        if (Vcycle == 2) then
                call plot('MGV02L4_D.txt',phi4,nx4,ny4,nz4)

        else if (Vcycle == 30) then
                call plot('MGV30L4_D.txt',phi4,nx4,ny4,nz4)
        else
        end if

        call restrict(res4,nx4,ny4,nz4,rhs5,nx5,ny5,nz5)

	call GS(5,nx5,ny5,nz5,4,phi5,rhs5,res5)

        if (Vcycle == 2) then
                call plot('MGV02L5_D.txt',phi5,nx5,ny5,nz5)

        else if (Vcycle == 30) then
                call plot('MGV30L5_D.txt',phi5,nx5,ny5,nz5)
        else
        end if

	call prolongate(phi54,nx4,ny4,nz4,phi5,nx5,ny5,nz5)
	phi4 = phi4 + phi54

	call GS(4,nx4,ny4,nz4,4,phi4,rhs4,res4)

        if (Vcycle == 2) then
                call plot('MGV02L4_U.txt',phi4,nx4,ny4,nz4)

        else if (Vcycle == 30) then
                call plot('MGV30L4_U.txt',phi4,nx4,ny4,nz4)
        else
        end if

        call prolongate(phi43,nx3,ny3,nz3,phi4,nx4,ny4,nz4)
        phi3 = phi3 + phi43

	call GS(3,nx3,ny3,nz3,6,phi3,rhs3,res3)

        if (Vcycle == 2) then
                call plot('MGV02L3_U.txt',phi3,nx3,ny3,nz3)

        else if (Vcycle == 30) then
                call plot('MGV30L3_U.txt',phi3,nx3,ny3,nz3)
        else
        end if

	call prolongate(phi32,nx2,ny2,nz2,phi3,nx3,ny3,nz3)
        phi2 = phi2 + phi32

        call GS(2,nx2,ny2,nz2,8,phi2,rhs2,res2)

        if (Vcycle == 2) then
                call plot('MGV02L2_U.txt',phi2,nx2,ny2,nz2)

        else if (Vcycle == 30) then
                call plot('MGV30L2_U.txt',phi2,nx2,ny2,nz2)
        else
        end if


	call prolongate(phi21,nx1,ny1,nz1,phi2,nx2,ny2,nz2)
        phi1 = phi1 + phi21

        call GS(1,nx1,ny1,nz1,10,phi1,rhs1,res1)

        if (Vcycle == 2) then
                call plot('MGV02L1_U.txt',phi1,nx1,ny1,nz1)

        else if (Vcycle == 30) then
                call plot('MGV30L1_U.txt',phi1,nx1,ny1,nz1)
        else
        end if

	Vcycle = Vcycle + 1
end do

close(200)

contains

subroutine plot(fname,phi,Nx,Ny,Nz)

! fname = file name
! phi is the variable which is to be plotted
! Nx,Ny,Nz = number of internal nodes for the plotted variable
! j = y-node at which the x-z contour is to be plotted
integer :: Nx,Ny,Nz,j
real,dimension(Nx+2,Ny+2,Nz+2) :: phi
character (len = 13) :: fname
	j = (Ny+2)/2
                open(unit = plot_var, file = fname)
                do k = 1,Nz+2,1
                write(plot_var,*) phi(:,j,k)
                end do
                close(plot_var)
plot_var = plot_var + 1

end subroutine plot


subroutine restrict(res_l1,Nx1,Ny1,Nz1,rhs_l2,Nx2,Ny2,Nz2)

! Nx1, Ny1, Nz1 are the number of internal nodes of res_l1 (fine)
! Nx2, Ny2, Nz2 are the number of internal nodes of rhs_l2 (coarse)

integer :: ic, jc, kc, i_f, jf, kf, Nx1, Ny1, Nz1, Nx2, Ny2, Nz2
real, dimension(Nx2+2,Ny2+2,Nz2+2), intent(out) :: rhs_l2
real, dimension(Nx1+2,Ny1+2,Nz1+2) :: res_l1

Do kc = 2,Nz2+1,1
	Do jc = 2,Ny2+1,1
		Do ic = 2,Nx2+1,1
			i_f = (ic-2) + ic 
			jf = (jc-2) + jc
			kf = (kc-2) + kc
			rhs_l2(ic,jc,kc) = 0.128*( res_l1(i_f,jf,kf) + res_l1(i_f+1,jf,kf) + res_l1(i_f+1,jf+1,kf) + res_l1(i_f,jf+1,kf) + &
				           res_l1(i_f,jf,kf+1) + res_l1(i_f+1,jf,kf+1) + res_l1(i_f+1,jf+1,kf+1) + res_l1(i_f,jf+1,kf+1) )
		end do
	end do
end do 

end subroutine restrict

subroutine prolongate(phi_21,Nx1,Ny1,Nz1,phi_2,Nx2,Ny2,Nz2)

! Nx1, Ny1, Nz1 are the number of internal nodes of phi_21 (fine)
! Nx2, Ny2, Nz2 are the number of internal nodes of phi_2 (coarse)

integer :: ic, jc, kc, i_f, jf, kf, Nx1, Ny1, Nz1, Nx2, Ny2, Nz2
real, dimension(Nx1+2,Ny1+2,Nz1+2), intent(out) :: phi_21
real, dimension(Nx2+2,Ny2+2,Nz2+2) :: phi_2

Do kf = 2,Nz1+2,2
	Do jf = 2,Ny1+2,2
		Do i_f = 2,Nx1+2,2
			ic = i_f/2 
			jc = jf/2
			kc = kf/2
			phi_21(i_f-1,jf-1,kf-1) = ( 27*phi_2(ic,jc,kc) + 9*phi_2(ic+1,jc,kc) + 9*phi_2(ic,jc,kc+1) + &
			9*phi_2(ic,jc+1,kc) + 3*phi_2(ic+1,jc+1,kc) + 3*phi_2(ic+1,jc,kc+1) + 3*phi_2(ic,jc+1,kc+1) + phi_2(ic+1,jc+1,kc+1) )/64
			
			phi_21(i_f,jf-1,kf-1) = ( 27*phi_2(ic+1,jc,kc) + 9*phi_2(ic+1,jc+1,kc) + 9*phi_2(ic,jc,kc) + &
			9*phi_2(ic+1,jc,kc+1) + 3*phi_2(ic+1,jc+1,kc+1) + 3*phi_2(ic,jc,kc+1) + 3*phi_2(ic,jc+1,kc) + phi_2(ic,jc+1,kc+1) )/64

			phi_21(i_f,jf,kf-1) = ( 27*phi_2(ic+1,jc+1,kc) + 9*phi_2(ic+1,jc,kc) + 9*phi_2(ic+1,jc+1,kc+1) + &
			9*phi_2(ic,jc+1,kc) + 3*phi_2(ic+1,jc,kc+1) + 3*phi_2(ic,jc+1,kc+1) + 3*phi_2(ic,jc,kc) + phi_2(ic,jc,kc+1) )/64

			phi_21(i_f-1,jf,kf-1) = ( 27*phi_2(ic,jc+1,kc) + 9*phi_2(ic,jc+1,kc+1) + 9*phi_2(ic,jc,kc) + &
			9*phi_2(ic+1,jc+1,kc) + 3*phi_2(ic+1,jc,kc) + 3*phi_2(ic+1,jc+1,kc+1) + 3*phi_2(ic,jc,kc+1) + phi_2(ic+1,jc,kc+1) )/64

			phi_21(i_f-1,jf-1,kf) = ( 27*phi_2(ic,jc,kc+1) + 9*phi_2(ic+1,jc,kc+1) + 9*phi_2(ic,jc+1,kc+1) + &
			9*phi_2(ic,jc,kc) + 3*phi_2(ic,jc+1,kc) + 3*phi_2(ic+1,jc,kc) + 3*phi_2(ic+1,jc+1,kc+1) + phi_2(ic+1,jc+1,kc) )/64

			phi_21(i_f,jf-1,kf) = ( 27*phi_2(ic+1,jc,kc+1) + 9*phi_2(ic+1,jc,kc) + 9*phi_2(ic,jc,kc+1) + &
			9*phi_2(ic+1,jc+1,kc+1) + 3*phi_2(ic+1,jc+1,kc) + 3*phi_2(ic,jc,kc) + 3*phi_2(ic,jc+1,kc+1) + phi_2(ic,jc+1,kc) )/64

			phi_21(i_f,jf,kf) = ( 27*phi_2(ic+1,jc+1,kc+1) + 9*phi_2(ic+1,jc+1,kc) + 9*phi_2(ic,jc+1,kc+1) + &
			9*phi_2(ic+1,jc,kc+1) + 3*phi_2(ic+1,jc,kc) + 3*phi_2(ic,jc,kc+1) + 3*phi_2(ic,jc+1,kc) + phi_2(ic,jc,kc) )/64

			phi_21(i_f-1,jf,kf) = ( 27*phi_2(ic,jc+1,kc+1) + 9*phi_2(ic,jc+1,kc) + 9*phi_2(ic,jc,kc+1) + &
			9*phi_2(ic+1,jc+1,kc+1) + 3*phi_2(ic+1,jc+1,kc) + 3*phi_2(ic+1,jc,kc+1) + 3*phi_2(ic,jc,kc) + phi_2(ic+1,jc,kc) )/64

		end do
	end do
end do 

end subroutine prolongate

subroutine GS(MGlevel,Nx,Ny,Nz,nits,phi_old,rhs,res)
! nits = No. of desired GS iterations at any MG level 
real :: del_x, del_y, del_z
real, dimension(Nx+2,Ny+2,Nz+2) :: phi, phi_diff, c, cw, ce, cs, cn, ca, cb, rhs
real, dimension(Nx+2,Ny+2,Nz+2), intent(out) :: phi_old, res
integer :: Nx, Ny, Nz, i, j, k, its, nits, MGlevel

! Number of internal nodes at any given MG level = Nx, Ny, Nz. Total nodes = Nx+2, Ny+2, Nz+2
! Grid size
del_x = Lx/(Nx)
del_y = Ly/(Ny)
del_z = Lz/(Nz)

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
do i = 2,Nx+1,1
        do j = 2,Ny+1,1
                do k = 2,Nz+1,1
                        c(i,j,k) = -2*( (1/(del_x*del_x))+(1/(del_y*del_y))+(1/(del_z*del_z)) )
                        cw(i,j,k) = 1/(del_x*del_x)
                        ce(i,j,k) = 1/(del_x*del_x)
                        cs(i,j,k) = 1/(del_y*del_y)
                        cn(i,j,k) = 1/(del_y*del_y)
                        ca(i,j,k) = 1/(del_z*del_z)
                        cb(i,j,k) = 1/(del_z*del_z)       
                end do
        end do
end do 

resd = 1
its = 0

! GS iterations
do while( its < nits )

do k = 1,Nz+2,1
        do j = 1,Ny+2,1
                do i = 1,Nx+2,1
                        
		if(MGlevel == 1) then
                        ! Boundary Conditions                
                        if (i == 1) then
                        phi(i,j,k) = 1000*exp(-50*(1 + (del_z*(0.5+k-2))**2)) - phi_old(i+1,j,k)
                        
                        else if (i == Nx+2) then
                        phi(i,j,k) = 2*( 100*(1-(del_z*(0.5+k-2))) + 500*exp(-50*((del_z*(0.5+k-2))**2))  ) - phi(i-1,j,k)
                        
                        else if (j == 1) then
                        phi(i,j,k) = phi_old(i,Ny+1,k)
                        
                        else if (j == Ny+2) then
                        phi(i,j,k) = phi(i,2,k)
                        
                        else if (k == 1) then
                        phi(i,j,k) = 2*( 100*(del_x*(0.5+i-2)) + 500*exp(-50*((1-(del_x*(0.5+i-2)))**2)) ) - phi_old(i,j,k+1)
                        
                        else if (k == Nz+2) then
                        phi(i,j,k) = 1000*exp(-50*(((1-(del_x*(0.5+i-2)))**2)+1)) - phi(i,j,k-1)

                        else
                        phi(i,j,k) = (rhs(i,j,k) - cw(i,j,k)*phi(i-1,j,k) - ce(i,j,k)*phi_old(i+1,j,k) - &
                        cs(i,j,k)*phi(i,j-1,k) - cn(i,j,k)*phi_old(i,j+1,k) - cb(i,j,k)*phi(i,j,k-1) - &
                        ca(i,j,k)*phi_old(i,j,k+1) )/c(i,j,k)                        

                        end if
		else

			if ((i == 1) .or. (i == Nx+2) .or. (j == 1) .or. (j == Ny+2) .or. (k == 1) .or. (k == Nz+2)) then

				phi(i,j,k) = 0
			else
				phi(i,j,k) = (rhs(i,j,k) - cw(i,j,k)*phi(i-1,j,k) - ce(i,j,k)*phi_old(i+1,j,k) - &
                        	cs(i,j,k)*phi(i,j-1,k) - cn(i,j,k)*phi_old(i,j+1,k) - cb(i,j,k)*phi(i,j,k-1) - &
                        	ca(i,j,k)*phi_old(i,j,k+1) )/c(i,j,k)
			end if
		end if                        
                end do
        end do
end do

phi_diff = phi - phi_old
resd = 0

! Calculating residual
do k = 2,Nz+1,1
        do j = 2,Ny+1,1
                do i = 2,Nx+1,1
                        resd = resd + (phi_diff(i,j,k))**2                        
                end do
        end do
end do

resd = sqrt(resd/(Nx*Ny*Nz))
phi_old = phi
its = its + 1

write(200,*) 'V cycle=',Vcycle, 'MG level=',MGlevel, 'GS its=',its, 'res=',resd


end do

do k = 2,Nz+1,1
        do j = 2,Ny+1,1
                do i = 2,Nx+1,1
                        res(i,j,k) = rhs(i,j,k) - cw(i,j,k)*phi(i-1,j,k) - ce(i,j,k)*phi(i+1,j,k) - &
                                cs(i,j,k)*phi(i,j-1,k) - cn(i,j,k)*phi(i,j+1,k) - cb(i,j,k)*phi(i,j,k-1) - &
                                ca(i,j,k)*phi(i,j,k+1) - c(i,j,k)*phi(i,j,k)
                end do
        end do
end do

!open(unit = 30, file = 'GS_N64.txt')
!do k = 1,Nz+2,1
!        write(30,*) phi(:,20,k)
!end do
!close(30)

end subroutine GS           

end program main
