#This example uses the banded function provided by the author to solve systems of equations involving a sparse matrix
#with a multidiagonal band.

import numpy as np
import banded as b
import pylab as py

 
bands=5

#specify number of unknowns, the amount of Voltage nodes we are calculating
N=6

A=np.zeros((bands,N))

v=np.zeros(N)

v[0]=5
v[1]=5



A[0,:]=-1
A[1,:]=-1
A[2,:]=4
A[2,0]=3
A[2,N-1]=3
A[3,:]=-1
A[4,:]=-1





V=np.round(b.banded(A,v,2,2),decimals=3);

i=np.linspace(1,6,6)


print('Voltage at each successive node starting from the first: ',V)

py.scatter(i,V)
py.title('Voltage at each node')

py.show()

#same solution using Gaussian Elimination and backsubsitution to Check

#Example 6.1

import numpy as np
A=np.array([[3,-1,-1,0,0,0],
            [-1,4,-1,-1,0,0],
           [-1,-1,4,-1,-1,0],
           [0,-1,-1,4,-1,-1],
           [0,0,-1,-1,4,-1],
           [0,0,0,-1,-1,3]],float)

v=np.array([5,5,0,0,0,0],float)

N=len(v)
for m in range(N):
    
    div=A[m,m]
    A[m,:]/=div
    v[m]/=div
    
    for i in range(m+1,N):
        mult=A[i,m]
        A[i,:]-=mult*A[m,:]
        v[i] -= mult*v[m]


x=np.empty(N,float)
for m in range(N-1,-1,-1):
    x[m]=v[m]
    for i in range(m+1,N):
        x[m] -=A[m,i]*x[i]
        


print(x)




##Now run for N=10,000 Voltage Nodes

bands=5

#specify number of unknowns, the amount of Voltage nodes we are calculating
N=10000

A=np.zeros((bands,N))

v=np.zeros(N)

v[0]=5
v[1]=5



A[0,:]=-1
A[1,:]=-1
A[2,:]=4
A[2,0]=3
A[2,N-1]=3
A[3,:]=-1
A[4,:]=-1


V=np.round(b.banded(A,v,2,2),decimals=3)

i=np.linspace(1,10000,10000)




py.scatter(i,V)
py.title('Voltage at each node')

py.show()
