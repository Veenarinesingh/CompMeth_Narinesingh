from math import sin
from numpy import arange,array,pi
from pylab import plot,xlabel,ylabel,show


g=9.81
l=.1

def f(r,t):
    theta=r[0]
    omega=r[1]
    ftheta=omega
    fomega=-(g/l)*sin(theta)
    return array([ftheta,fomega],float)

a = 0.0
b = 10.0
N = 1000

h = (b-a)/N
tpoints = arange(a,b,h)
thetapoints = []
omegapoints=[]

r=array([3.124,0],float)

for t in tpoints:
    thetapoints.append(r[0])
    omegapoints.append(r[1])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3=h*f(r+.5*k2,t+.5*h)
    k4=h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6



plot(tpoints,thetapoints)
xlabel("t")
ylabel("x(t)")
show()
