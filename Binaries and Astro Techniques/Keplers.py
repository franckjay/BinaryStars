import math

def Keplers(e,MA):
    
    E = MA
    
    ecc=float(e)
    
    M = MA
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
      
    return E

def trueAnomaly(e,E):
    a=(1+e)
    b=(1-e)
    c=math.tan(E/2)        
    theta = 2*math.atan(math.sqrt(a/b)*c)
    
    return theta
    
    
