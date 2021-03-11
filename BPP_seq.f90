! AUTH. ABHIJIT BAISHYA
! LAST MODIFIED 09.03.2021

program pec_Seq_decay
implicit none

double precision :: m(2), m12, z(2), AN(2), R0, rad(2), hbar, amu, R1, R2, aa, bb, d2r, r2d
double precision :: Vmin, kmin, k, S, l, Ex, Eth, Ecm, P, T

integer :: i, j, No
double precision time, start, finish

    call cpu_time(start)


d2r=0.017453292
r2d=1./d2r


Ex = 7.654       !Hoyle State excitation energy of the recoil 12C
Eth = 7.274      !Three alpha breakup threshold
Ecm = Ex - Eth   !Exc. energy in C.M.

hbar = 197.326 !Mev-fm
amu = 931.494 !MeV


!Mass , indices for 3 different alpha
m(1) = 8.0
m(2) = 4.0

m12 = m(1)*m(2)/(m(1)+m(2)) !Reduced Mass of the System

!A of Nuclei , indices for 3 different alpha
AN(1) = 8.0
AN(2) = 4.0


!Charge of Nuclei , indices for 3 different alpha
z(1) = 4.0
z(2) = 2.0

open(10,file='BPP_seq.txt')

do i=1,16

R0 = 0.8+i*0.05 ! in fm

l = 0.0  ! Angular Momentum, e.g. for Hoyle State is 0+, l is 0


! fragment radii
rad(1) = R0*(AN(1)**(1.0/3.0))
rad(2) = R0*(AN(2)**(1.0/3.0))




P = 0.0


!print *,'Excitation Energy:'
!print *,Ex


! ######### Penetrability Calculation Starts Here   ########


R1 = rad(1) + rad(2)
!print *,R1

aa = R1
bb = 100.0 !! As solve uses BISECTION METHOD, you need to provide to initial points, 1st one is rho0, 2nd one is > rho0

call solve(aa,bb,R2,m12,z,l,Ecm)
!write(*,*)R2

No = 10000  ! Integration Points for the subroutine simpsint
call simpsint(R1, R2, No, m12, z, l, Ecm, S)

kmin = sqrt(2.0*m12*amu*f(R1,m12,z,l,Ecm))/hbar
k = sqrt(2.0*m12*amu*Ecm)/hbar

S = S*(sqrt(2.0*m12*amu)/hbar)
!print *,S


!P = exp(-2.0*S) ! Alpha Decay Probability
!P = (k/(2.0*R1))*exp(-2.0*S) ! Alpha Decay Rate
P = kmin/(1.0+exp(2.0*S))  ! R-matrix Probability

T = (3.0*10.0**(-23.0))*(2.0*R1/sqrt(2.0*Ecm/m12))*exp(2.0*S) ! Lifetime

print *,'Probability for Tunneling:'
print *,P

print *,'Decay Time Constant:'
print *,T, "sec"
write(10,*)P
enddo
! ######### Penetrability Calculation Ends Here   ########


! Code Run Time
    call cpu_time(finish)
    time = finish - start

print *, 'CPU run time'
print *, time, 'sec'


close(10)




!!!! ##### Subroutines and Functions #####
contains   !All the required functions and subroutines




!Function for root solve "(V - Exc.)=0"
function f(x,m,z,l,Ex)
implicit none

double precision :: f, x, m, z(2), l, Ex, hbar

hbar = 197.326

f = (1.44)*(z(1)*z(2)/x) + ((hbar**2.0)*l*(l+1))/(2.0*m*amu*x**2.0) - Ex
!f = x**2.0 - 1.0

return
end function 






!!Simpson's 1/3 rule integration
subroutine simpsint(a, b, N, m, z, l, Ex, res)
implicit none 
double precision :: a, b, res, h, m, z(2), l, Ex
integer :: i, N
h = (b-a)/N
res = 0.0
do i=1,N-1
res = res + (h/6.0)*(sqrt(f(a+(i-1)*h,m,z,l,Ex)) + 4.0*sqrt(f(a+(i-0.5)*h,m,z,l,Ex)) + sqrt(f(a+i*h,m,z,l,Ex)))
enddo

return
end subroutine










!!Subroutine for root solving
subroutine solve(a,b,x,m,z,l,Ex)
double precision :: a, b, x, temp, m, z(2), l, Ex
integer :: i, j, N

do i=1,100000
if (f(a,m,z,l,Ex)*f(b,m,z,l,Ex) .lt. 0.) then

	x = (a+b)/2.0
        temp = x
	if (f(a,m,z,l,Ex)*f(x,m,z,l,Ex) .lt. 0.) then

		a = a
		b = x
	elseif (f(b,m,z,l,Ex)*f(x,m,z,l,Ex) .lt. 0.) then

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
if (abs(f(x,m,z,l,Ex)) .lt. 1.0e-7) exit
enddo

return
end subroutine






end program

