import math
import Keplers
import Periods2
import matplotlib.pyplot as plt

e=.350
P=43
To=1992.3
delta = P/36.0
w=math.radians(290)
inc= math.radians(107)
a = 0.21
Omega=math.radians(164)
t = []
phi = []
MA=[]
E=[]
r=[]
ValX=[]
ValY=[]
theta=[]


for i in range(36):
    a = To + (i*delta)
    t.append(a)

for i in t:
    a=(Periods2.orbitalPhase(43,1992.3,float(i)))
    phi.append(a)


for i in phi:
    b = i*2*math.pi
    MA.append(b)

for i in MA:
    c = Keplers.Keplers(e,i)
    E.append(c)

for i in E:
    d = a*(1-(e*math.cos(i)))    
    theta.append(Keplers.trueAnomaly(e,i))
    r.append(d)

cnt=0
for i in r:
    x=i*((math.cos(Omega))*(math.cos(theta[cnt]+w))-((math.sin(Omega))*(math.sin(theta[cnt]+w))*(math.cos(inc))))
    y=i*((math.sin(Omega))*(math.cos(theta[cnt]+w))+((math.cos(Omega))*(math.sin(theta[cnt]+w))*(math.cos(inc))))
    ValX.append(x)
    ValY.append(y)
    cnt+=1


    
plt.plot([ValY],[ValX], 'ro')
plt.show()


    



    




    
