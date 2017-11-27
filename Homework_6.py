from vpython import *
from numpy import arange


ball = sphere(pos=vector(0,0,0), radius=0.5, color=color.cyan,make_trail=True)

r=vector(0,0,0)

for n in arange(0,1e6):
    rate(1e5)
    ball.pos=r
    if r.x==50 and r.y==50:
        a=random()
        if a<1/2:
            r.x=r.x-1
        if a>1/2:
            r.y=r.y-1
    elif r.x==50 and r.y==-50:
        a=random()
        if a<1/2:
            r.x=r.x-1
        if a>1/2:
            r.y=r.y+1
    elif r.x==-50 and r.y==-50:
        a=random()
        if a<1/2:
            r.x=r.x+1
        if a>1/2:
            r.y=r.y+1
    elif r.x==-50 and r.y==50:
        a=random()
        if a<1/2:
            r.x=r.x+1
        if a>1/2:
            r.y=r.y-1
    elif r.x==50:
        a=random()
        if a<1/3:
            r.x=r.x-1
        if 1/3<a<2/3:
            r.y=r.y+1 
        if 2/3<a:
            r.y=r.y-1
    elif r.x==-50:
        a=random()
        if a<1/3:
            r.x=r.x+1
        if 1/3<a<2/3:
            r.y=r.y+1 
        if 2/3<a:
            r.y=r.y-1
    elif r.y==50:
        a=random()
        if a<1/3:
            r.x=r.x+1
        if 1/3<a<2/3:
            r.x=r.x-1 
        if 2/3<a:
            r.y=r.y-1
    elif r.y==-50:
        a=random()
        if a<1/3:
            r.x=r.x+1
        if 1/3<a<2/3:
            r.x=r.x-1 
        if 2/3<a:
            r.y=r.y+1
            
    else:
        b=random()
        if b<.25:
            r.x=r.x+1
        elif .25<b<.50:
            r.x=r.x-1
        elif .50<b<.75:
            r.y=r.y+1
        elif .75<b<1:
            r.y=r.y-1
        
    

   
