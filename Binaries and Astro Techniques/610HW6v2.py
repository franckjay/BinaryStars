#########################################################################
## Jay Franck
## AST 610 Homework 6 #2
## Newton-Raphson Iteration for Kepler's Equation
#########################################################################

import math
import Keplers
import matplotlib.pyplot as plt
import pylab

def syntheticVelocity(gamma,K,e,w,theta):
    V= gamma + K*((e*math.cos(w))+ (math.cos(theta+w)))
    return V


eccSpu1 = .12
gamma = 30
K = 35
phi=0
w = 0
phiList=[]
vSpu1= []
obsVel= []
obsPhase= []
obsSpu1= open("outputSpu1.txt", 'r')
for i in obsSpu1:
    i=i.split()
    obsVel.append(i[0])
    obsVel.append(i[0])
    obsPhase.append(i[1])
    obsPhase.append(float(i[1])+1)

for i in range(401):
    phiList.append(phi)
    ##Loop through each phi with each ecc and w, finding the V, Store V
    if phi >= 1:
        M = (phi-1)*(2)*(math.pi)
    else:
        M = (phi)*(2)*(math.pi)
    E=Keplers.Keplers(eccSpu1,M)
   
    theta=Keplers.trueAnomaly(eccSpu1,E)
    

    vSpu1.append(syntheticVelocity(gamma,K,eccSpu1,w,theta))

    phi+= 0.005


plt.xlabel('Phase')
plt.ylabel('Synthetic Velocities (km/s)')
plt.title('Velocity Curves for Spu 1')
plt.plot(phiList,vSpu1,'r-')
plt.plot(obsPhase,obsVel,'ro')

plt.show()


