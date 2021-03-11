# AUTH. ABHIJIT BAISHYA
# LAST MODIFIED 09.03.2021
# Program for Penetrability for sequential decay having bent arm structure

import numpy as np
from matplotlib import pyplot as plt
#from sympy import symbols, solve
#from random import random
import scipy.integrate as integrate



class pec:

    d2r=0.017453292
    r2d=1./d2r

    def __init__(self, AN, R0, Z, Ex, th):
        self.AN = AN
        self.Z = Z
        self.Ex = Ex
        self.R0 = R0
        self.th = th

    def distance(self, r):
        self.r = r
        rho = np.zeros(3)
        rad = np.zeros(3)
        rad[0] = self.R0*(self.AN[0]**(1.0/3.0))
        rad[1] = self.R0*(self.AN[1]**(1.0/3.0))
        rad[2] = self.R0*(self.AN[2]**(1.0/3.0))


        #calculate the distances between the particles
        r12sq = (rad[0] + rad[1])**2.0
        if self.r == 0.0 :
            r23sq = (rad[1] + rad[2])**2.0
        else:
            r23sq = self.r**2.0
        r13sq = r12sq + r23sq - 2.0*np.sqrt(r12sq*r23sq)*np.cos(self.th*pec.d2r) #-th[2])

        #the hyperradius of the system
        #Mt = np.sum(m)

        #rhosq = (1/(sm*Mt))*(m[0]*m[1]*r12sq + m[0]*m[2]*r13sq + m[1]*m[2]*r23sq)
        #rho = np.sqrt(rhosq)
        rho[0] = np.sqrt(r12sq)
        rho[1] = np.sqrt(r23sq)
        rho[2] = np.sqrt(r13sq)
        return rho

    def Vcoulcf(self,r):
        self.r = r
        rho = pec.distance(self,r)
	    #V = 1.44*( Z[0]*Z[1]/rho[0] + Z[1]*Z[2]/rho[1] +  Z[2]*Z[0]/rho[2] ) - Ex
        V = 1.44*( self.Z[1]*self.Z[2]/rho[1] +  self.Z[2]*self.Z[0]/rho[2] ) - self.Ex
        return V

    
    def eqsolve(self,aa,bb):
        self.aa = aa
        self.bb = bb 
        x = (self.aa+self.bb)/2.
        while (abs(pec.Vcoulcf(self,self.aa) - pec.Vcoulcf(self,self.bb)) > 1.0e-06):
            if (pec.Vcoulcf(self,self.aa)*pec.Vcoulcf(self,self.bb) <= 0.) :
                x = (self.aa+self.bb)/2.
                temp = x
                if (pec.Vcoulcf(self,self.aa)*pec.Vcoulcf(self,x) <= 0.) :
                	self.aa = self.aa
                	self.bb = x
                else:
                	self.aa = x
                	self.bb = self.bb
            else:
            	self.bb = 2.*self.bb

        return (self.aa+self.bb)/2.
   


hbar = 197.326 #Mev-fm
amu = 931.494 #MeV

    
Exc = 7.65       #Hoyle State excitation energy of the recoil 12C
#Exc(2) = 10.03      #1st Excited state of Hoyle State for the recoil 12C
Eth = 7.36      #Breakup threshold

Ex = Exc - Eth

Pl = []
ang_arm = []
f = open("pec_seq_bent_arm.txt", 'w')

for i in range(36):
    #Mass , indices for 3 different alpha
    m = [4.0, 8.0]
    #A of Nuclei , indices for 3 different alpha
    AN = [4.0, 4.0, 4.0]
    #Charge of Nuclei , indices for 3 different alpha
    Z = [2.0, 2.0, 2.0]

    R0 = 1.05 # in fm

    #sm = 3.0/8.0 #Hyperspherical normalizing mass
    m12 = (m[0]*m[1])/(m[0]+m[1]) #Reduced Mass

    ang_arm.append((i+1)*5.0)

    alp_Be = pec(AN, R0, Z, Ex, (i+1)*5.0)
    
    rho0 = alp_Be.distance(0.0)
    #print(rho0)
    #print(alp_Be.Vcoulcf(0.))
    r0 = rho0[1]
    r1 = alp_Be.eqsolve(r0,100.)
    #print(r1,alp_Be.Vcoulcf(r1))


    N = 500
    h = (r1-r0)/N
    x = np.zeros(N+1)

    for j in range(0,N):
        x[j] = r0 + j*h

    integrand = np.zeros(N+1)
    for k in range(0,N):
        integrand[k] = np.sqrt(alp_Be.Vcoulcf(x[k]))
    
    S = integrate.trapz(integrand, x) 
    
    kmin = np.sqrt(2*m12*amu*(alp_Be.Vcoulcf(r0)))/hbar
    k = np.sqrt(2*m12*amu*Ex)/hbar

    S = S*(np.sqrt(2*m12*amu)/hbar)
    #print(S)

    P = kmin*np.exp(-2.0*S)  #R-matrix Probability
    print('Probability for Tunneling:', P)
    f.write('{}'.format((i+1)*5.0) + "  ")
    f.write('{}'.format(P) + "\n")
    Pl.append(P)

f.close()

plt.plot(ang_arm,Pl,label='Bent Arm Angle Dependence of Penetrability')
plt.xlabel('Angle of Bent Arm in degrees', size=10)
plt.ylabel('Penetrability', size=10)
#plt.show()
