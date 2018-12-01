import math


c=299792.5

def radVel(lamM,lamLab):

    V = c*(lamM-lamLab)/lamLab

    return V

def meanV(vList):
    sumVel= 0
    for i in vList:
        i = float(i)
        sumVel+=i
    mean= sumVel/len(vList)

    return mean

def residuals(Vm,mean):
    deviation = Vm - mean
    return deviation
    
def stdev(deviations,vList):
    num = 0
    N = len(vList)
    ##FIX STDEV FOR INPUT OF RESIDUALS
    for i in deviations:
        OC2 = math.pow(i,2)
        num += OC2
    a = num/(N-1)
    SD = math.sqrt(a)
    return SD

lamMea=[6347.178,6371.473,6456.382,6516.345,6562.944,6678.307,6402.409]
lamLab=[6347.095,6371.355,6456.277,6516.053,6562.808,6678.152,6402.250]
Species=['Si II','Si II','Fe II + H20?','Fe II','Halpha','He I','Ne I ??']

vList= []
deviations=[]
mean=0
standardDev=0

for i in range(len(lamLab)):
    vList.append(radVel(lamMea[i],lamLab[i]))
mean = meanV(vList)

for i in vList:
    deviations.append(residuals(i,mean))

standardDev = stdev(deviations,vList)
helioRV = mean+0.14-6.53
cnt = 1
for i in vList:
    if i > (mean +(3*standardDev)):
        print "Reject measurement: ",cnt," due to Chauvenet's criterion!"
    elif i < (mean -(3*standardDev)):
        print "Reject measurement: ",cnt," due to Chauvenet's criterion!"

print "The Mean is: ",mean," km/s"
print "The Standard Deviation is: ",standardDev," km/s"
print "The Heliocentric Radial Velocity is ",helioRV," km/s"
print "____________________________________________________________________"
print "Species||Lab Wave||Meas Wave||Velocity||Deviation"
print "____________________________________________________________________"
for i in range(len(lamLab)):
    print Species[i],lamLab[i],lamMea[i],vList[i],deviations[i]




