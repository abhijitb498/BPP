program Coulomb_Penetribility
!  use omp_lib

implicit none

double precision :: E1, Ex, Exc(2), Eth, Ecm1, Ecm2, Ecm3, E3, E4, theta3, &
 theta4, phi3, phi4, P3, P4, d2r, r2d, a , b , d , r , s, t, u
double precision :: Ealp1, Ealp2, Ealp3, E_Be, v1(3), v_1(3), v2(3), v_2(3), &
v3(3), v_3(3), v_Be(3), v__Be(3), v_C(3) 
double precision :: theta_1, phi_1, theta_2, phi_2, theta_3, phi_3, theta, phi, psi, jai
double precision :: th_phi1(2), th_phi2(2), th_phi3(2)
double precision :: e12, e23, e31, xdalitz, ydalitz
double precision :: v11,v22,v33, factor, theta_phi3(2)

double precision :: m(3), z(3), AN(3), sm, kl, kl1(2), rad(3), th(3), v(3), rho, rho0, rho1, ss(3), s0(3), aa, bb
double precision :: Vcc, Vmin, kmin, k, S1, Pl, Pb, P1(2), P(2)
!real, allocatable :: P1(:), P(:)
integer :: No, Ns, Np
!For Inputs
integer :: i, j, N, N_hit_seq, N_hit_DDL, N_hit_DDE
double precision :: theta_min, theta_max, phi_min, phi_max
double precision :: rho_range

double precision wtime, start, finish
    call cpu_time(start)


d2r=0.017453292
r2d=1./d2r


m(1) = 4.0
m(2) = 4.0
m(3) = 4.0

AN(1) = 4.0
AN(2) = 4.0
AN(3) = 4.0

z(1) = 2.0
z(2) = 2.0
z(3) = 2.0

sm = 1.0 !Hyperspherical normalizing mass

! hypermomentum quantum number
Kl1(1) = 0.0
kl1(2) = 2.0


! fragment radii
rad(1) = 1.05*(AN(1)**(1.0/3.0))
rad(2) = 1.05*(AN(2)**(1.0/3.0))
rad(3) = 1.05*(AN(2)**(1.0/3.0))


!write(*,*)'ENTER THE NUMBER OF EVENTS:'
!read(*,*)N
!N = 5000


Exc(1) = 7.654       !Hoyle State excitation energy of the recoil 12C
Exc(2) = 10.03      ! 1st Excited state of Hoyle State for the recoil 12C
Eth = 7.274      ! Three alpha breakup threshold

!Calling Random Seed to initialize the Random Number Generator
call randomEG(r)

!!!!##### START OF MAIN LOOP #######
!P1 = 0.0

! !$OMP PARALLEL DO
!do j=1,2
Ns = 0
Ex = Exc(1)
kl = kl1(1)
!print *,'Excitation Energy:'
!print *,Ex

! !$OMP PARALLEL DO
!do i=1,N




!########### DECAY IN DDphi MODE ###########




Ecm1 = Ex - Eth

!Distributing energy between 3 alphas'
!call random_number(r)
!call random_number(s)
v11 = 1./3
v22 = 1./3
v33 = 1.0/3.0 ! - v22 - v11

!if (v11+v22 .ge. v33 .and. v22+v33 .ge. v11 .and. v33+v11 .ge. v22) then   ! Momentum in CM zero means three vectors should follow Triangle rule
!Ns = Ns + 1
!In CMF emission angles
call random_number(r)
theta_1 = 90.0   !Polar Angle of emitted alpha-1
phi_1 = 360.*r     !Azimuthal Angle of emitted alpha-1

th(1) = acos((v33**2.0 - v11**2.0 - v22**2.0)/(2.*v11*v22))  !Angle bet. alp-1 & alp-2

theta_2 = 90.0
phi_2 = phi_1 + th(1)*r2d

th(2) = acos((v11**2.0 - v22**2.0 - v33**2.0)/(2.*v33*v22))  !Angle bet. alp-2 & alp-3

theta_3 = 90.0
phi_3 = phi_2 + th(2)*r2d

th(3) = acos((v22**2.0 - v11**2.0 - v33**2.0)/(2.*v33*v11))  !Angle bet. alp-1 & alp-3

v_1 = v11*pol_cart(theta_1,phi_1)
v_2 = v22*pol_cart(theta_2,phi_2)
v_3 = v33*pol_cart(theta_3,phi_3)


factor = sqrt(Ecm1/(energy(m(1), v_1) + energy(m(2), v_2) + energy(m(3), v_3)))

v_1 = factor*v_1
v_2 = factor*v_2
v_3 = factor*v_3

!write(*,*)(v_1 + v_2 + v_3)

v11 = sqrt(v_1(1)**2.0 + v_1(2)**2.0 + v_1(3)**2.0)
v22 = sqrt(v_2(1)**2.0 + v_2(2)**2.0 + v_2(3)**2.0)
v33 = sqrt(v_3(1)**2.0 + v_3(2)**2.0 + v_3(3)**2.0)

!write(*,*)(energy(4., v_1) + energy(4., v_2) + energy(4., v_3))
v(1) = v11
v(2) = v22
v(3) = v33


!Potential for DDE Config
Np = 200 !No. of points for the plot
rho_range = 140 !Range of plot
call hyperradius2(m, sm, th, v, rho, ss) !Find out the scaling factors s(i)
!write(*,*)ss
call hyperradius1(m, sm, th, rad, rho0, s0) !Finds out the initial hyperradius rho0
!write(*,*)rho0

rho = rho0
 open(10,file='potential_DDE.txt')
 open(11,file='potential_hoyle.txt')
 open(12,file='potential_2+.txt')

rho = rho0-2.5
do i=1,Np
Vcc = Vcoulcf(rho,sm,z,kl,ss)
!write(10,*)rho0, Vcc, 0.38, rho
write(10,*)rho, Vcc
write(11,*)rho, 0.38
write(12,*)rho, 2.756
!rho0 = rho0 + 0.1
!rho = rho + (rho_range-rho0)/Np
rho = rho + 1.0
end do

 close(10)
 close(11)
 close(12)

!Potential for DDL Config 
call hyperradius_DDL(m, sm, rad, rho0, s0)
rho = 14.37
!print*,Vcoulcf(rho,sm,z,kl,s0)
open(13,file='potential_DDL.txt')

rho = rho0-2.5
do i=1,Np
Vcc = Vcoulcf(rho,sm,z,kl,s0)
write(13,*)rho, Vcc
rho = rho + 1.0
!rho = rho + (rho_range-rho0)/Np
end do

 close(13)

 

!call execute_command_line('gnuplot -p potential.plt')

!!!! ##### Subroutines and Functions #####
contains   !All the required functions and subroutines



!! Polar to Cartesian Conversion
function pol_cart(theta,phi)
double precision :: pol_cart(3)
double precision :: theta, phi, d2r, r2d 

d2r=0.017453292
r2d=1./d2r

theta = theta*d2r
phi = phi*d2r

pol_cart(1) = sin(theta)*cos(phi)
pol_cart(2) = sin(theta)*sin(phi)
pol_cart(3) = cos(theta)
return 
end function






!!Find Theta and Phi from vector Components
function pol_angles(vec)
double precision :: pol_angles(2)
double precision :: vec(3)
double precision :: theta, phi, d2r, r2d, r 

d2r=0.017453292
r2d=1./d2r
r = sqrt(vec(1)**2. + vec(2)**2. + vec(3)**2.)

pol_angles(1) = acos(vec(3)/r)*r2d
pol_angles(2) = atan(vec(2)/vec(1))*r2d
return
end function





!!Find Energy of a particle given its velocity vector
function energy(mass, vec)
double precision :: energy
double precision :: vec(3), mass

energy = 0.5*mass*(vec(1)**2. + vec(2)**2. + vec(3)**2.)
return
end function





!!Dalitz Plot Coordinate Transformation
function dalitz(e12, e23, e31)
double precision :: dalitz(2)
double precision :: e12, e23, e31

dalitz(1) = sqrt(3.)*(e12 - e23)
dalitz(2) = (2.*e31 - e12 - e23)
return
end function







!!Coulomb + Centrifugal potential
function Vcoulcf(rho,sm,z,k,s)
implicit none
double precision :: Vcoulcf
double precision :: rho,sm,z(3),k,s(3)

Vcoulcf = (1.44/rho)*(z(1)*z(2)/s(1) + z(2)*z(3)/s(2) + z(1)*z(3)/s(3)) + &
(197.326**2.0)*(k+1.5)*(k+2.5)/(2.0*sm*931.494*rho**2.0) - 132 / (1 + exp((-1.7 + rho*s(1)) / 0.669)) &
- 132 / (1 + exp((-1.7 + rho*s(2)) / 0.669)) - 132 / (1 + exp((-1.7 + rho*s(3)) / 0.669)) !in MeV when rho is in fm and sm in MeV
return
end function





!Function for root solve
function f(x,sm,z,k,s,Ex)
implicit none

double precision :: f,x,sm,z(3),k,s(3), Ex

f = (1.44/x)*(z(1)*z(2)/s(1) + z(2)*z(3)/s(2) + z(1)*z(3)/s(3)) + &
((197.326**2.0)*(k+1.5)*(k+2.5))/(2.0*sm*931.494*x**2.0) - Ex
!f = x**2.0 - 1.0

return
end function 






!!!! Subroutine for calculating initial hyperradius
subroutine hyperradius1(m, sm, th, D, rho, s)
implicit none
double precision :: m(3), sm, th(3), D(3)
double precision :: rho, s(3), Mt, rhosq, r12sq, r13sq ,r23sq, d2r, r2d

d2r=0.017453292
r2d=1./d2r


!calculate the distances between the particles
!r12sq = D(1)**2.0 + D(2)**2.0 - 2*D(1)*D(2)*cos(th(1))
!r13sq = D(1)**2.0 + D(3)**2.0 - 2*D(1)*D(3)*cos(th(2))
!r23sq = D(2)**2.0 + D(3)**2.0 - 2*D(2)*D(3)*cos(th(3))
r12sq = (D(1) + D(2))**2.0
r23sq = (D(2) + D(3))**2.0
r13sq = r12sq + r23sq - 2.0*sqrt(r12sq*r23sq)*cos(th(3))
!r13sq = (D(1) + D(3))**2.0

!the hyperradius of the system
Mt = m(1) + m(2) + m(3)

rhosq = (1/(sm*Mt))*(m(1)*m(2)*r12sq + m(1)*m(3)*r13sq + m(2)*m(3)*r23sq)
rho = sqrt(rhosq)

!calculate the scaling parameters
s(1) = sqrt(r12sq/rhosq)
s(2) = sqrt(r23sq/rhosq)
s(3) = sqrt(r13sq/rhosq)
return
end subroutine





!! Subroutine for calculating scaling factors
subroutine hyperradius2(m, sm, th, D, rho, s)
implicit none
double precision :: m(3), sm, th(3), D(3)
double precision :: rho, s(3), Mt, rhosq, r12sq, r13sq ,r23sq, d2r, r2d

d2r=0.017453292
r2d=1./d2r

!calculate the distances between the particles
r12sq = D(1)**2.0 + D(2)**2.0 - 2*D(1)*D(2)*cos(th(1))
r23sq = D(2)**2.0 + D(3)**2.0 - 2*D(2)*D(3)*cos(th(2))
r13sq = D(1)**2.0 + D(3)**2.0 - 2*D(1)*D(3)*cos(th(3))

!the hyperradius of the system
Mt = m(1) + m(2) + m(3)

rhosq = (1/(sm*Mt))*(m(1)*m(2)*r12sq + m(1)*m(3)*r13sq + m(2)*m(3)*r23sq)
rho = sqrt(rhosq)

!calculate the scaling parameters
s(1) = sqrt(r12sq/rhosq)
s(2) = sqrt(r23sq/rhosq)
s(3) = sqrt(r13sq/rhosq)
return
end subroutine







!! Subroutine for calculating hyperradius-2
subroutine hyperradius_DDL(m, sm, r, rho, s)
implicit none
double precision :: m(3), sm, r(3)
double precision :: rho, s(3), Mt, rhosq, r12sq, r13sq ,r23sq

!calculate the distances between the particles
r12sq = (4.0*r(1))**2.0
r13sq = (2.0*r(2))**2.0
r23sq = (2.0*r(3))**2.0

!the hyperradius of the system
Mt = m(1) + m(2) + m(3)

rhosq = (1/(sm*Mt))*(m(1)*m(2)*r12sq + m(1)*m(3)*r13sq + m(2)*m(3)*r23sq)
rho = sqrt(rhosq)

!calculate the scaling parameters
s(1) = sqrt(r12sq/rhosq)
s(2) = sqrt(r13sq/rhosq)
s(3) = sqrt(r23sq/rhosq)
return
end subroutine




!!Simpson's 1/3 rule integration
subroutine simpsint(a, b, N, sm, z, k, s, Ex, res)
implicit none 
double precision :: a,b,res,h, sm, z(3), k, s(3), Ex
integer :: i, N
h = (b-a)/N
res = 0.0
!  !$OMP PARALLEL DO
do i=1,N-1
if (f(a,sm,z,k,s,Ex) .lt. 0. .or. f(a+(i-0.5)*h,sm,z,k,s,Ex) .lt. 0.0 .or. f(a+(i)*h,sm,z,k,s,Ex) .lt. 0.0) exit !then
res = res + (h/6.0)*(sqrt(f(a,sm,z,k,s,Ex)) + 4.0*sqrt(f(a+(i-0.5)*h,sm,z,k,s,Ex)) + sqrt(f(a+i*h,sm,z,k,s,Ex)))
!print *,(h/6.0)*(sqrt(f(a,sm,z,k,s,Ex)) + 4.0*sqrt(f(a+(i-0.5)*h,sm,z,k,s,Ex)) + sqrt(f(a+i*h,sm,z,k,s,Ex)))
!print *,f(a+i*h,sm,z,k,s,Ex)
!endif
enddo
! !$OMP END PARALLEL DO

return
end subroutine







!!Subroutine for root solving
subroutine solve(a,b,x,sm,z,k,s,Ex)
double precision :: a,b,x,temp,sm,z(3),k,s(3), Ex
integer :: i,j,N
! !$OMP PARALLEL DO

do i=1,100000
if (f(a,sm,z,k,s,Ex)*f(b,sm,z,k,s,Ex) .lt. 0.) then

	x = (a+b)/2.0
        temp = x
	if (f(a,sm,z,k,s,Ex)*f(x,sm,z,k,s,Ex) .lt. 0.) then

		a = a
		b = x
	elseif (f(b,sm,z,k,s,Ex)*f(x,sm,z,k,s,Ex) .lt. 0.) then

		a = x
		b = b
	endif
!	write(*,*)x,temp
!	if (abs(f(x,sm,z,k,s,Ex)-f(temp,sm,z,k,s,Ex)) .lt. 1.0e-4 .or. f(x,sm,z,k,s,Ex) .eq. f(temp,sm,z,k,s,Ex)) exit
!write(*,*)a,b
else
	 b = 2*b
endif
!write(*,*)x
if (abs(f(x,sm,z,k,s,Ex)) .lt. 1.0e-7) exit
enddo
! !$OMP END PARALLEL DO

return
end subroutine








!Random Seed Generator
subroutine randomEG(r)
implicit none

integer :: k, i, n=10
double precision :: r
integer, dimension(8) :: values 
! Declare an assumed shape, dynamic array
integer, dimension(:), allocatable :: seed

! gfortran subroutine to return date and time information 
! from the real time system clock. Works down to milliseconds 
! and stores the eight return values in array values.
call date_and_time(VALUES=values)
! restart the state of the pseudorandom number generator
! k = minimum size of seed (12 on my system)
call random_seed(size=k)
! allocate memory to seed
allocate(seed(k))

! assign information in values to seed
seed(:) = values(:)
! seed the random number generator
call random_seed(put=seed)

!do i=1,n
!    call random_number(r)
!    print *, r
!end do
return
end subroutine randomEG

end program
