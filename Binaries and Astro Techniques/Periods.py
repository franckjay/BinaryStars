#########################################################################
##Jay Franck
##AST 610 Homework 2
#########################################################################
import md5
import random
import math

def orbitalPhase(P,To,t):
##This subroutine will recieve a time from the file HW02data and will find
##the orbital phase and the cycle number.
    P=float(P)
    To=float(To)
    t=float(t)
    ##This converts all values to double floating precision values
    print "p: ",P
    print "T: ",T
    print "t: ",t    
    phi =((t-To))/P
    print phi
    if phi < 0:
        phi += 1
    while phi > 1:
        ##This loop changes To to the starting position of each new orbit
        n=2
        To = To + (P*n)
        phi = (abs(t-To)/P)
        n +=1
    ##I calculates the cycle number for each period and converts it to an int
    I = int(abs((t-6912.034)/P))-1
    ##This corrects the fact that the cycle does not start until one full
    ##period from To.
    #period = open('period.txt', 'a')  
    #c = "____________________________________________________________________"
    #a=str("Orbital Phase: %.7f Cycle Number: %i Time: %f" %(phi,I,t))
    #period.write(a+"\n")
    #period.write(c+"\n")
    #period.close()
   
    return phi

##This is the main body of the program where subroutines are called##
##The data file is called, and then the for loop iterates through each time
##calling the subroutine individually
times = open('HW02data.txt', 'r')
P = 1.1164072
To = 6912.034
for i in times:
    orbitalPhase(P,To,i)
    
    
    
