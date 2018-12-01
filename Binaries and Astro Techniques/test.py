import math
import numpy
import mpfit

def myfunct(p,fjac=None,x=None,y=None,err=None):
    model = p[0] +p[1]*x+ p[2]*(x**2) +p[3]*numpy.sqrt(x)+p[4]*numpy.log(x)
    status=0
    return ([status, (y-model)])

err=1
x = numpy.arange(101, dtype=float)
p = [5, 2, 500, 1,2000]
y = (p[0]+ p[1]*x + p[2]*(x**2) +p[3]*numpy.sqrt(x)+p[4]*numpy.log(x) )
fa = {'x':x, 'y':y, 'err':err}
m = mpfit(myfunct(p,x,y), p, functkw=fa)
print 'status = ', m.status
if (m.status <= 0):
    print 'error message = ', m.errmsg
print 'parameters = ', m.params
   

