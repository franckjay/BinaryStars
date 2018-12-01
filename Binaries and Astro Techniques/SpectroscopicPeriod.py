#########################################################################
##Jay Franck
##AST 610 Homework 4: Finding Periods from minimal data
#########################################################################

import Periods2
import matplotlib.pyplot as plt
import math


times=[]
radVel=[]
trialPeriods=[]
theta=[]
thetaBest=1000000000.0
periodBest=0

dataFile = open('1Spu.txt', 'r')
for i in dataFile:
    a = i.split()
    times.append(a[0])
    radVel.append(a[1])

##This is the number of observations
N=len(times)
To = times[0]

## den(denominator) is the sum of the squared difference of the mean radVel(M) 
## from each individual radVel.
a = 0
for i in radVel:
    a+=float(i)
M=a/len(radVel)
den=0
for i in radVel:
    difference=(float(i)-M)
    den+=math.pow(difference,2)

   
##Creates a list of trial periods
for i in range(10,11):
    for j in range(1,1000):
	    k=float(j)
	    a = float(i)+float(k/1000)
	    trialPeriods.append(a)


cnt= 0
while cnt < len(trialPeriods):
    cnt2=0
    tempVel=radVel
    phase = []
    P =trialPeriods[cnt]
    
    ##Iterates through each time, calculating the phase, append it to a list
    while cnt2 < len(times):
        a = Periods2.orbitalPhase(P,To,times[cnt2])
        phase.append(a)
        cnt2+=1

    ##Quick and dirty parallel sorter: Use Python's zip function to iterate
    ##through two sets of lists at once, then sort the phases using the
    ##built-in Python adaptive_merge_sort algorithm. Then use a 'map' function
    ##to with arbitrary map function 'lambda' that correlates the value in list1
    ##to list2, thus sorting the array.
    data=zip(phase,tempVel)
    data.sort()
    phase,tempVel=map(lambda t: list(t), zip(*data))
    
    ##Discriminant Function from Lafler and Kinman
    numerator=0
    cnt3=0

    while cnt3 < 28:
        cnt3=int(cnt3)
        diff = float(tempVel[cnt3])-float(tempVel[cnt3+1])
        numerator += math.pow(diff,2)
        cnt3+=1

    tempTheta = numerator/den
    theta.append(tempTheta)
    if tempTheta < thetaBest:
        thetaBest=tempTheta
        periodBest=P
    cnt +=1    
    
out = open('outputSpu1.txt', 'a')

    
##Taking the PeriodBest as an input, the program gathers the phases once again
##to plot a radial velocity versus phase diagram
phasesBest =[]
for i in range(len(times)):
    bestPhase = Periods2.orbitalPhase(10.33538,To,times[i])
    phasesBest.append(bestPhase)
    w = str(radVel[i])+" "+str(bestPhase)+"\n"
    out.write(w)
plt.plot([phasesBest],[radVel], 'ro')
plt.ylabel('radVel')
plt.xlabel('Phase')
plt.show()
out.close()
##This plots the disciminating function versus the Period in days
#plt.plot([trialPeriods],[theta], 'ro')
#plt.ylabel('Theta')
#plt.xlabel('P in days')
#plt.show()



