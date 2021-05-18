! AUTH. ABHIJIT BAISHYA
! LAST MODIFIED 09.03.2021

program Penetribility_DDphi
!  use omp_lib
 use ieee_arithmetic
implicit none

double precision :: E1, Ex, Exc(2), Eth, Ecm1, Ecm2, Ecm3, E3, E4, theta3, &
 theta4, phi3, phi4, P3, P4, d2r, r2d, a , b , d , r , s, t, u
double precision :: Ealp1, Ealp2, Ealp3, E_Be, v1(3), v_1(3), v2(3), v_2(3), &
v3(3), v_3(3), v_Be(3), v__Be(3), v_C(3) 
double precision :: theta_1, phi_1, theta_2, phi_2, theta_3, phi_3, theta, phi, psi, jai
double precision :: th_phi1(2), th_phi2(2), th_phi3(2)
double precision :: e12, e23, e31, xdalitz, ydalitz
double precision :: v11,v22,v33, factor, theta_phi3(2)

double precision :: m(3), z(3), AN(3), R0, sm, kl, kl1(2), rad(3), hbar, amu, th(3), v(3), rho, rho0, rho1, ss(3), s0(3), aa, bb
double precision :: Vmin, kmin, k, S1, Pl, Pb, P1(2), P(2), P_seq
double precision :: max_val
character*20 :: files(2)
integer :: No, Ns
!For Inputs
integer :: i, j, ii, N, N_hit_seq, N_hit_DDL, N_hit_DDE
double precision :: theta_min, theta_max, phi_min, phi_max

double precision time, start, finish
    call cpu_time(start)


d2r=0.017453292
r2d=1./d2r

Exc(1) = 7.654       !Hoyle State excitation energy of the recoil 12C
Exc(2) = 10.03      ! 1st Excited state of Hoyle State for the recoil 12C
Eth = 7.274      ! Three alpha breakup threshold

hbar = 197.326 !Mev-fm
amu = 931.494 !MeV

!Mass , indices for 3 different alpha
m(1) = 4.0
m(2) = 4.0
m(3) = 4.0

!A of Nuclei , indices for 3 different alpha
AN(1) = 4.0
AN(2) = 4.0
AN(3) = 4.0

!Charge of Nuclei , indices for 3 different alpha
z(1) = 2.0
z(2) = 2.0
z(3) = 2.0

! hypermomentum quantum number
Kl1(1) = 0.0 !For Hoyle State
kl1(2) = 2.0 !For Hoyle 1st Excited

sm = 1.0 !Hyperspherical normalizing mass

files = (/"potential_DDphi.txt","potential_DDph2.txt"/)
!open(10,file='BPP_DDphi.txt')
!open(11,file='BPP_seq.txt')

write(*,*)'ENTER THE NUMBER OF EVENTS:'
!read(*,*)N
N = 1000
!Calling Random Seed to initialize the Random Number Generator
call randomEG(r)


!do ii=1,16

!R0 = 0.8+ii*0.05 ! in fm
R0 = 1.05 !in fm

!read(11,*)P_seq

! fragment radii
rad(1) = R0*(AN(1)**(1.0/3.0))
rad(2) = R0*(AN(2)**(1.0/3.0))
rad(3) = R0*(AN(3)**(1.0/3.0))


!!!!##### START OF MAIN LOOP #######
P1 = 0.0

do j=1,2   ! Starts Loop for two different Excitation Energies, e.g. 7.654 MeV & 10.03 MeV
Ns = 0
Ex = Exc(j)
kl = kl1(j)
print *,'Excitation Energy:'
print *,Ex

i = 1
do while (i .le. N)  ! Starts Loop for different Phase Space Configurations

!########### DECAY IN DDphi MODE ###########

Ecm1 = Ex - Eth

!Distributing energy between 3 alphas'
call random_number(r)
call random_number(s)
v11 = 0.5*r
v22 = (1-v11)*s
v33 = 1.0 - v22 - v11

v = rand_shuffle((/v11,v22,v33/))
v11 = v(1)
v22 = v(2) 
v33 = v(3)

if (v11+v22 .ge. v33 .and. v22+v33 .ge. v11 .and. v33+v11 .ge. v22) then   ! Momentum in CM zero means three vectors should follow Triangle rule
!Ns = Ns + 1
i = i+1
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

!!########### BARRIER PENETRIBILITY CALCULATION STARTS HERE ############

call hyperradius2(m, sm, th, v, rho, ss) !Find out the scaling factors s(i)
!write(*,*)ss
call hyperradius1(m, sm, th, rad, rho0, s0) !Finds out the initial hyperradius rho0
!print*, "Barrier at rho0", rho0, "is: ", f(rho0,sm,z,kl,ss,Ecm1) + Ecm1

!Vmin = Vcoulcf(rho0,sm,z,kl,ss)

!if (minval(th)*r2d .lt. 26.0) then
!print*, minval(th)*r2d
!open(13,file=files(j))
!rho = rho0
!do i=1,100
!!print*, f(rho,sm,z,kl,ss,Ecm1) + Ecm1
!write(13,*)rho, f(rho,sm,z,kl,ss,Ecm1) + Ecm1
!rho = rho + 0.5
!end do
!close(13)
!endif

call get_max(max_val)
if (Ecm1 .lt. max_val) then
aa = rho0
bb = rho0 !! As solve uses BISECTION METHOD, you need to provide to initial points, 1st one is rho0, 2nd one is > rho0
call solve(aa,bb,rho0,sm,z,kl,ss,Ecm1)
!print*,"rho0 is: ", rho0
!print*, "Barrier at rho0", rho0, "is: ", f(rho0,sm,z,kl,ss,Ecm1) + Ecm1


aa = rho0*1.000001
bb = aa
call solve(aa,bb,rho1,sm,z,kl,ss,Ecm1)
!write(*,*)rho1

No = int(rho1/0.05)  ! Integration Points for the subroutine simpsint

rho0 = rho0*1.000001
call simpsint(rho0, rho1, No, sm, z, kl, ss, Ecm1, S1)
S1 = S1*(sqrt(2*sm*amu)/hbar)

!print*,"rho0, rho1 & S1", rho0, rho1, S1
!print*, "Barrier at rho0 = ", rho0, "is", f(rho0,sm,z,kl,ss,Ecm1) + Ecm1
	
Pl = 1.0/(1.0+exp(2.0*S1)) ! Alpha Decay Probability
!Pl = (k/(2.0*rho0))*exp(-2.0*S1) ! Alpha Decay Rate
!Pl = k*exp(-2.0*S1)
!print*, Pl
Ns = Ns + 1
else
!print*, max_val, th*r2d
Pl = 1.0
Ns = Ns + 1
endif

kmin = sqrt(2*sm*amu*(Vcoulcf(rho0,sm,z,kl,ss) - Ecm1))/hbar
!print*, "kmin is: ", kmin
k = sqrt(2*sm*amu*Ecm1)/hbar
Pl = kmin*Pl  ! R-matrix Probability
!Pl = k*Pl
!print*, Pl
P1(j) = P1(j) + Pl
!endif
! ##### BARRIER PENETRIBILITY CALCULATION ENDS ###########


endif

!End Of 3rd Do loop for Phase Space Configurations
enddo

print *,'Probability for Tunneling:'
P(j) = P1(j)/Ns
print *,P(j)


enddo ! End of 2st do loop

print *,'Ratio of the two Probablities for R0: ', R0
print *, P(1)/P(2)

!write(10,*)R0, P, P(1)/P(2), P(1)/P_seq

!enddo ! End of 1st looop

!close(10)
!close(11)

call cpu_time(finish)
time = finish - start

print *, 'CPU run time'
print *, time, 'sec'










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
if (vec(1) .gt. 0.0 .and. vec(2) .gt. 0.0) then
	pol_angles(2) = atan(vec(2)/vec(1))*r2d + 0.0
!phi = acos(vec(1)/(r*sin(acos(vec(3)/r))))*r2d
elseif (vec(1) .lt. 0.0 .and. vec(2) .gt. 0.0) then
	pol_angles(2) = atan(vec(2)/vec(1))*r2d + 180.0
elseif (vec(1) .lt. 0.0 .and. vec(2) .lt. 0.0) then
	pol_angles(2) = atan(vec(2)/vec(1))*r2d + 180.0
else
	pol_angles(2) = atan(vec(2)/vec(1))*r2d + 360.0
endif
!pol_angles(2) = atan(vec(2)/vec(1))*r2d
return
end function





!!Find Energy of a particle given its velocity vector
function energy(mass, vec)
double precision :: energy
double precision :: vec(3), mass

energy = 0.5*mass*(vec(1)**2. + vec(2)**2. + vec(3)**2.)
return
end function






!!Coulomb + Centrifugal potential
function Vcoulcf(rho,sm,z,k,s)
implicit none
double precision :: Vcoulcf
double precision :: rho,sm,z(3),k,s(3), hbar

hbar = 197.326

Vcoulcf = (1.44/rho)*(z(1)*z(2)/s(1) + z(2)*z(3)/s(2) + z(1)*z(3)/s(3)) + &
(hbar**2.0)*(k+1.5)*(k+2.5)/(2.0*sm*amu*rho**2.0) !in MeV when rho is in fm and sm in MeV
return
end function





!Function for root solve Vcoulcf - Excitation Energy
function f(x,sm,z,k,s,Ex)
implicit none

double precision :: f,x,sm,z(3),k,s(3), Ex, hbar
double precision :: V_d, r_d, a_d

hbar = 197.326
V_d = -132.0
r_d = 1.7
a_d = 0.669

!print*, "s is: ", s
f = (1.44/x)*(z(1)*z(2)/s(1) + z(2)*z(3)/s(2) + z(1)*z(3)/s(3)) + &
((hbar**2.0)*(k+1.5)*(k+2.5))/(2.0*sm*amu*x**2.0) + V_d / (1 + exp((-r_d + x*s(1)) / a_d)) &
+ V_d / (1 + exp((-r_d + x*s(2)) / a_d)) + V_d / (1 + exp((-r_d + x*s(3)) / a_d)) - Ex
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






!!Simpson's 1/3 rule integration
subroutine simpsint(a, b, N, sm, z, k, s, Ex, res)
implicit none 
double precision :: a,b,res,temp,h, sm, z(3), k, s(3), Ex
integer :: i, N
h = (b-a)/N
res = 0.0
!  !$OMP PARALLEL DO
do i=1,N-1
!print*,(f(a+(i-1)*h,sm,z,k,s,Ex))
temp = (h/6.0)*(sqrt(f(a+(i-1)*h,sm,z,k,s,Ex)) + 4.0*sqrt(f(a+(i-0.5)*h,sm,z,k,s,Ex)) + sqrt(f(a+i*h,sm,z,k,s,Ex)))

if (ieee_is_nan(temp) .eqv. .false.) then
res = res + (h/6.0)*(sqrt(f(a+(i-1)*h,sm,z,k,s,Ex)) + 4.0*sqrt(f(a+(i-0.5)*h,sm,z,k,s,Ex)) + sqrt(f(a+i*h,sm,z,k,s,Ex)))
endif

enddo

return
end subroutine





!!Subroutine for root solving
subroutine solve(a,b,x,sm,z,k,s,Ex)
double precision :: a,b,x,temp,sm,z(3),k,s(3), Ex
integer :: i,j,N
! !$OMP PARALLEL DO

do !i=1,10000
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
if (abs(f(x,sm,z,k,s,Ex)) .lt. 1.0e-07) exit
enddo
! !$OMP END PARALLEL DO

return
end subroutine


subroutine get_max(max)
implicit none
double precision, intent(out) :: max
double precision :: arr(400), rho
integer :: i

rho = 2.0

do i=1,400
arr(i) = f(rho,sm,z,kl,ss,Ecm1)
rho = rho + 0.5
enddo

max = maxval(arr)

return
end subroutine

!Randomly Shuffle an array of 3 elements
function rand_shuffle(arr)
implicit none
double precision :: rand_shuffle(3), r, a, b, c
double precision, intent(in) :: arr(3)
integer :: i

a = arr(1)
b = arr(2)
c = arr(3)

!call randomEG(r)
call random_number(r)

i = int(r*6)
!print*,i

select case(i)
    case(0)
        rand_shuffle = (/a,b,c/)
    case(1)
        rand_shuffle = (/a,c,b/)
    case(2)
        rand_shuffle = (/b,a,c/)
    case(3)
        rand_shuffle = (/c,a,b/)
    case(4)
        rand_shuffle = (/b,c,a/)
    case(5)
        rand_shuffle = (/c,b,a/)
end select

return
end function rand_shuffle



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
