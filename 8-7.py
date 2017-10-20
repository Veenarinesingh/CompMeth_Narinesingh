import numpy as np
import pylab as py

#Constants
R=.08
rho=1.22
C=.47
theta_in_degrees=30
theta=theta_in_degrees*np.pi/180
g=9.8
for m in range(1,4):

    if m==1:
        d='C1'
    if m==2:
        d='C2'
    if m==3:
        d='C3'


    D=np.pi*(R**2)*rho*C/(2*m)

    def f(r,t):
        x=r[0]
        xdot=r[1]
        y=r[2]
        ydot=r[3]
        fx=xdot
        fy=ydot
        fxdot=-D*xdot*np.sqrt(xdot**2+ydot**2)
        fydot=-g-D*ydot*np.sqrt(xdot**2+ydot**2)
        return np.array([fx,fxdot,fy,fydot])

    a=0.0
    b=10
    N=1000

    h=(b-a)/N

    tpoints = np.arange(a,b,h)
    xpoints = []
    xdotpoints=[]
    ypoints=[]
    ydotpoints=[]

    r=np.array([0,100*np.cos(theta),0,100*np.sin(theta)])

    for t in tpoints:
        xpoints.append(r[0])
        xdotpoints.append(r[1])
        ypoints.append(r[2])
        ydotpoints.append(r[3])
        k1 = h*f(r,t)
        k2 = h*f(r+0.5*k1,t+0.5*h)
        k3=h*f(r+.5*k2,t+.5*h)
        k4=h*f(r+k3,t+h)
        r += (k1+2*k2+2*k3+k4)/6

    py.plot(xpoints,ypoints,d,label=m)


py.axis([min(xpoints),max(xpoints),0,max(ypoints)+10])
py.xlabel("Horizontal Distance (meters)")
py.ylabel("Height (meters)")
py.legend(title='Mass of canonball in kg')
py.show()

#Heavier objects have less air resistance, and thus travel further
