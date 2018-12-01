#########################################################################
## Jay Franck
## AST 680 Homework 4
## Newton-Raphson Iteration for Kepler's Equation
#########################################################################
import math

##Introduce variables and add the list of Mean Anomalies
pi=math.pi
ecc=0.666
M = [(pi/4),(pi/2),(3*pi/4)]

##Iterates through each value of M and creates the header file
for i in M:
    i = float(i)
    print ""
    print "Newton-Raphson Iteration of Kepler's Equation"
    print "Paramters are M: ",i,"e: ",ecc
    print ""
    print "%4s %10s %10s %10s %10s %10s"%("i",'Ei','dEi','Ei+1','num','den')
    print "_______________________________________________________________"
    N=1
    temp=0
    E=i
    ##Starts the NR Iteration with the set parameters
    while N<10:
        f = -i+E-(ecc*(math.sin(E)))
        df = 1-(ecc*(math.cos(E)))
        tmp = E - (f/df)            
        dE = abs(E-tmp)
        print "%4i %10f %10f %10f %10f %10f"%(N,E,dE,tmp,f,df)
        E=tmp
        N+=1
        if dE < 0.000001:
            break
        
    
