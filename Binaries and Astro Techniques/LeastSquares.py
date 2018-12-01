import math
import matplotlib.pyplot as plt

times =[]
timesFile = open('minTimes.txt', 'r')
for i in timesFile:
    i = i.strip('\t\n')
    times.append(i)

To = 37968.3438
E = []
N=len(times)
##The maximum possible period must be at least the minimum
##distance between any two minimum times. 
minimum=100000000000
for i in range(N-1):
    m = float(times[int(i+1)])-float(times[int(i)])
    if m < minimum:
        minimum = m
minimum=2*minimum
minimum = math.ceil(minimum)

##Accounting for any possible errors in observation times,
##spurios periods, etc. the possible range of periods will be
##double the minimum value.
trialPeriods = []
#Narrow down the period
for i in range(2,int(minimum)):
    for j in range(1,20000):
	    k=float(j)
	    a = float(i)+float(k/20000)
	    trialPeriods.append(a)

nTP=len(trialPeriods)
minChi=100000000
periodBest= 0
##Performs the Chi squared test through the list of trial periods
##saving the period with the minimum value of Chi as a good approx
##of the actual period
for i in range(nTP):
    i = int(i)
    Chi = 0
    for k in times:
        k = float(k)
        E=int((float(k)-To)/trialPeriods[i])
        Chi += math.pow((k-(To+(E*trialPeriods[i]))),2)
    if Chi < minChi:            
        minChi = Chi
        periodBest = trialPeriods[i]

##Find Ba from the periodBest
Ei=[]
for i in times:
    j=int((float(i)-To)/periodBest)
    Ei.append(j)
xiyi=0
xi = 0
squareXi=0
yi=0
for i in range(N):
    i=int(i)
    E = float(Ei[i])
    t = float(times[i])
    xiyi += (E*t)
    xi += E
    yi += t
    squareXi += math.pow(E,2)
Ba = ((N*xiyi)-(xi*yi))/((N*squareXi)-math.pow(xi,2))
##Find Aa from Ba
nextCycle=[]
for i in times:
    j=int((float(i)-To)/Ba)
    nextCycle.append(j)
xiyi=0
xi = 0
squareXi=0
yi=0
for i in range(N):
    i=int(i)
    E = float(nextCycle[i])
    t = float(times[i])
    xiyi += (E*t)
    xi += E
    yi += t
    squareXi += math.pow(E,2)
Aa = ((squareXi*yi)-(xiyi*xi))/((N*squareXi)-math.pow(xi,2))
       
Ebest = []
OC =[]
aChi = 0
for k in times:
    k = float(k)
    E=int((float(k)-Aa)/Ba)
     Ebest.append(E)
    oc = (k-(Aa+(E*Ba)))
    OC.append(oc)
    aChi += math.pow((k-(Aa+(E*Ba))),2)

##Error Calculations
##Based on best values
meanX =0
for i in Ebest:
    i=float(i)
    meanX +=i
meanX = meanX/N

SSxx = 0
for i in Ebest:
    i=float(i)
    SSxx += math.pow((i-meanX),2)
    
meanY= 0
for i in times:
    i=float(i)
    meanY += i
meanY =meanY/N

SSyy = 0
for i in times:
    i=float(i)
    SSyy += math.pow(i,2)
SSyy = SSyy -(N*math.pow(meanY,2))

SSxy = 0
for i in range(N):
    i = int(i)
    SSxy += float(times[i])*float(Ebest[i])
SSxy = SSxy -(N*meanX*meanY)

s= math.sqrt((SSyy - ((math.pow(SSxy,2))/SSxx))/(N-2))

sigA = math.sqrt((1/N)+(math.pow(meanX,2)/SSxx))
sigB = s/math.sqrt(SSxx)

plt.plot([Ebest],[OC], 'ro')
plt.ylabel('O-C')
plt.xlabel('Cycle Number E')
plt.title('O-C Diagram of a binary using Linear Regression')
plt.show()


