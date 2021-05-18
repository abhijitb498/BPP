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


!Potential for DDE Config
Np = 400 !No. of points for the plot
ss = (/0.29708700047081421, 0.56138540127795555, 0.58871533476141358/)
open(10,file='potential_DDphi.txt')
open(11,file='potential_DDphi_2.txt')

rho = 5.0
do i=1,Np

write(10,*)rho, Vcoulcf(rho,sm,z,kl1(1),ss)
write(11,*)rho, Vcoulcf(rho,sm,z,kl1(2),ss)

rho = rho + 0.5
end do

close(10)
close(11)


contains


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


end program
