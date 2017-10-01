#Exercise 5.3

import numpy as np

import pylab as py


#integrating the function


#define the function

def f(t):
    return np.exp(-t**2)


#show a plot of what the function looks like
r=np.linspace(0,3,31)
g=f(r)

#print out some helpful details for the user

print('f(t)=exp(-t'u"\u00B2"')')
py.plot(r,g)
py.xlabel('t')
py.ylabel('f(t)')
py.show()



print('E(x) is the integral of f(t) from 0 to x')

#Use trapezoidal rule to evaluate integral, a is starting point, b=x is ending point,h is width of slices, N is number of slices
#create an x vector for the range we are interested in and an empty y vector to be populated by our integral value.

xvector=np.linspace(0,3,31)
yvector=[]


N=1000
a=0


for n in range(0,31):
    
    x=xvector[n]

    h=(x-a)/N

    s=.5*f(a)+.5*f(x)

    for k in range(1,N):
        s+=f(a+k*h)
        
    print('Integral of f(t) from 0 to',x,'is',h*s)
    
    yvector=yvector+[h*s]

    
#plot the results
py.plot(xvector,yvector)
py.xlabel('x')
py.ylabel('E(x)')

py.show()





