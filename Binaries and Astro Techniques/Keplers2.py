#########################################################################
##Jay Franck
##AST 610 Homework 3 Problem #1
#########################################################################
import math
from decimal import *
getcontext().prec = 7

def Keplers(e):
    meanAnomaly = [0,45,90,135,180,225,270,315]
    e = [0.00,0.15,0.5,0.85]
    E = math.radians(math.pi)
    OutputE= []
    OutputdE=[]
    OutputN=[]
    for a in e:
0000000000000000        ecc=float(a)
        for i in meanAnomaly:
            M = math.radians(i)
            N=1
            tmp = 0
            while N < 25:
                f = -M+E-(ecc*(math.sin(E)))
                df = 1-(ecc*(math.cos(E)))
                tmp = E - (f/df)            
                dE = abs(E-tmp)
                E=tmp
                N+=1
                if dE < 0.00001:
                    break
            OutputE.append(E)
            OutputdE.append(dE)
            OutputN.append(N)
    
        
    
    return OutputE,OutputdE,OutputN

def trueAnomaly(e,E):
    thetaOut = []
    cnt = 0

    for i in e:
        a=(1+e[cnt])
        b=(1-e[cnt])
        c=math.tan(E[cnt]/2)
        
        theta = 2*math.atan(math.sqrt(a/b)*c)
        thetaOut.append(theta)
        cnt+=1
    
    return thetaOut


e =[0.00,0.15,0.5,0.85]
meanAnomaly = [0,45,90,135,180,225,270,315]
N = []
Edeg = []
Erad= []
dEdeg=[]
dErad=[]
cnt = 0
eOut=[]
MAout=[]
MArad=[]
rOut = []
TArad=[]
TAdeg=[]

Keplers(e)

for i in Keplers(e)[0]:
    Erad.append(i)
    Edeg.append(math.degrees(i))
for i in Keplers(e)[1]:
    dErad.append(i)
    dEdeg.append(math.degrees(i))
for i in Keplers(e)[2]:
    N.append(i)
for i in e:
    for j in meanAnomaly:
        eOut.append(i)
        MAout.append(j)
        MArad.append(math.radians(j))
        r=1-(i*math.cos(Erad[cnt]))
        rOut.append(r)
        cnt+=1
for i in trueAnomaly(eOut, Erad):
    TArad.append(i)
    TAdeg.append(math.degrees(i))

    
    
period = open('Anomalies.txt', 'a')
Header="%s %10s %10s %10s %10s %10s %4s %12s %12s %5s"%('e','M(deg)','M(rad)','E(deg)','E(rad)','de(rad)','N','theta(deg)','theta(rad)','r')
period.write(Header+'\n')
period.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+'\n')

for i in range(32):
    data="%.2f %10.6f %10.6f %10.6f %10.6f %10.6f %i %12.6f %10.6f %10.6f"%(eOut[i],MAout[i],MArad[i],Edeg[i],Erad[i],dErad[i],N[i],TAdeg[i],TArad[i],rOut[i])
    period.write(data+'\n')
    
period.close()
        

