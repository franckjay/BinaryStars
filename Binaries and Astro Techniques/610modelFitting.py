#########################################################################
## Jay Franck
## AST 610 Homework 7 
## Model Fitting
#########################################################################

import math
import Periods2
import Keplers
#import numpy
#import mpfit
#import matplotlib.pyplot as plt


def syntheticVelocity(gamma,K,e,w,theta):
    V=gamma+K*((e*math.cos(w))+(math.cos(theta+w)))
    return V


##Uses the values from the Spumoni data and its calculated orbital phase
##from Homework 6 and adds them to a list
obsVel= []
times= []
obsSpu1= open("1Spu.txt", 'r')
for i in obsSpu1:
    i=i.split()
    obsVel.append(i[1])
    times.append(i[0])
for i in range(len(times)):
    obsVel.append(obsVel[i])
    
##
bestP=0
bestTo=0
beste=0
bestw=0
bestK=0
bestGam=0
##
minChi=100000000
P = [10.3365]
To = [45785.665]
for z in range(5):
    con = .05
    To.append(float(times[0])+con)
    To.append(float(times[0])-con)  
    con +=.05
e=[0.008]
w=[0.002]
K=[33.835]
gamma=[33.40]



for a in P:    
    for b in To:
        ##Calculates the orbital phase for each point
        obsPhase=[]
        for i in times:
            phase = Periods2.orbitalPhase(a,b,i)
            obsPhase.append(phase)
        for i in range(len(obsPhase)):
                       phase = obsPhase[i]+1
                       obsPhase.append(phase)              
        for c in e:
            theta=[]
            for i in obsPhase:
                ##Calculates Mean Anomaly, Eccentric Anomaly, then True Anomaly
                i = float(i)
                if i >= 1:
                    M = (i-1)*(2)*(math.pi)
                else:
                    M = (i)*(2)*(math.pi)
                E=Keplers.Keplers(c,M)
                t=Keplers.trueAnomaly(c,E)
                theta.append(t)
            for d in w:
                for f in K:
                    for g in gamma:
                        Chi=0
                        ##Chi Squared test
                        for j in range(len(theta)):
                            v=syntheticVelocity(g,f,c,d,theta[j])
                            ov=float(obsVel[j])
                            Chi += math.pow((ov-v),2)
                    if Chi < minChi:
                        minChi=Chi
                        bestP=a
                        bestTo=b
                        beste=c
                        bestw=d
                        bestK=f
                        bestGam=g


a = open('radVelout.txt', 'a')
for i in range(len(obsPhase)):
    a.write(str(obsPhase[i])+" "+str(obsVel[i])+'\n')

a.close()

phi=0
phiList=[]
vSpu1=[]
print minChi
for i in range(401):
    phiList.append(phi)
    if phi >= 1:
        M = (phi-1)*(2)*(math.pi)
    else:
        M = (phi)*(2)*(math.pi)
    E=Keplers.Keplers(beste,M)
   ##Creates line of best fit
    theta=Keplers.trueAnomaly(beste,E)
    vSpu1.append(syntheticVelocity(bestGam,bestK,beste,bestw,theta))

    phi+= 0.005

#plt.xlabel('Phase')
#plt.ylabel('Synthetic Velocities (km/s)')
#$plt.title('Velocity Curves for Spu 1')
#plt.plot(obsPhase,obsVel,'ro')
#plt.plot(phiList,vSpu1,'b--')

#plt.show()        

                        




