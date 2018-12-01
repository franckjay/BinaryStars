#########################################################################
## Jay Franck
## AST 610 Homework 5
## Newton-Raphson Iteration for Kepler's Equation
#########################################################################

import math
import Keplers
import matplotlib.pyplot as plt


def syntheticVelocity(gamma,K,e,w,theta):
    V= gamma + K*((e*math.cos(w))+ (math.cos(theta+w)))
    return V


ecc1 = 0.7
ecc2 = 0.02
gamma = 25
K = 50
phi=0
phiList=[]
##List of omegas to iterate through, with empty lists of synthetic RadVels
w1 = [0,90,180,270]
V0=[]
V90=[]
V180=[]
V270=[]


w2 = [30,120,210,300]
V30=[]
V120=[]
V210=[]
V300=[]

for i in range(401):
    phiList.append(phi)
    ##Loop through each phi with each ecc and w, finding the V, Store V
    if phi >= 1:
        M = (phi-1)*(2)*(math.pi)
    else:
        M = (phi)*(2)*(math.pi)
    E1=Keplers.Keplers(ecc1,M)
    E2=Keplers.Keplers(ecc2,M)
    t1=Keplers.trueAnomaly(ecc1,E1)
    t2=Keplers.trueAnomaly(ecc2,E2)

    V0.append(syntheticVelocity(gamma,K,ecc1,w1[0],t1))
    V90.append(syntheticVelocity(gamma,K,ecc1,w1[1],t1))
    V180.append(syntheticVelocity(gamma,K,ecc1,w1[2],t1))
    V270.append(syntheticVelocity(gamma,K,ecc1,w1[3],t1))

    V30.append(syntheticVelocity(gamma,K,ecc2,w2[0],t2))
    V120.append(syntheticVelocity(gamma,K,ecc2,w2[1],t2))
    V210.append(syntheticVelocity(gamma,K,ecc2,w2[2],t2))
    V300.append(syntheticVelocity(gamma,K,ecc2,w2[3],t2))
    
    phi+= 0.005


plt.xlabel('Phase')
plt.ylabel('Synthetic Velocities (km/s)')
plt.title('e=0.02 Velocity Curves for the family of w=30,120,210,300 degrees')
plt.plot(phiList,V30,'g--')    
plt.plot(phiList,V120,'r.')
plt.plot(phiList,V210,'b,')
plt.plot(phiList,V300,'b^')

plt.show()
