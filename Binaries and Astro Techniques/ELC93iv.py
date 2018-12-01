##Homework 9, part 3iii
import math


def X(w):
    w=math.radians(w)
    a = ((0.3*math.cos(w))/(math.sqrt(1-math.pow(0.3,2))))
    X= math.pi + 2*math.atan(a)
    return X



w = [0,30,60,90,120,150,180,210,240,270,300,330]
a = open('ELC3iv.txt', 'a')
for i in w:
    i = float(i)
    c= X(i)
    b = c - math.sin(c)
    b = b/(2*math.pi)
    a.write(str(i) + " " + str(b)+ "\n")

a.close()
    

