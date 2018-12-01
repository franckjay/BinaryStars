import math




def R(R):

    
    N=1
    tmp = 0
    crit = 9.47E-27
    matter= 0.27*crit
    rel = 8.24E-5*crit
    
    c1 = (8/3)*math.pi*6.67E-11
    c2 = (0.044*crit*2.997E8*6.65E-29)/1.67E-27
    
    while N < 100:
        f = c1*((math.pow(R,3)*matter)+(math.pow(R,2)*rel))-math.pow(c2,2)
        df = (c1*3*math.pow(R,2)*matter)+(c1*2*R*rel)
        tmp = R - (f/df)            
        dR = abs(R-tmp)
        R=tmp
        N+=1
        if dR < 0.00001:
            break
        print N
        

    return R

print R(.01)
