#Veeshan Narinesingh, Computational Methods in Physics Fall '17
#Professor Ari Maller

#Final Project
#Numerical Solution to the Barotropic Vorticity Equation, Initial Condition of
#single negative Gaussian depression, Periodic Boundary in x and y

#Note:This code has been adapted from Peter Lynch's MatLab code which can be
#accessed here https://maths.ucd.ie/met/msc/fezzik/MatLab/matlab.html


#    d              g                       d
#   -- (Del^2-F)w + - J(w,Del^2(w)) + beta* -- w = 0.
#   dt              f                      dx
#
#     With a mean zonal flow term
#     Ubar*d((Del^2w)/dx added.

from numpy import linspace,pi,sin,cos,sqrt,meshgrid,size,shape,zeros,transpose,fix,mean,array,real,asarray,float32,exp
from random import seed,random
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import time as Time
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from numpy.fft import fft2,ifft2

#Initial Conditions Specifications for a pseudo-real 500 mb flow
daylen=input('How long of a forecast? (integers only please!)')
DAYLEN=int(daylen)         # Forecast length in days.
NX = 50             # Set spatial resolution
NY = 50
DELTA_t = 1/24      # Timestep in hours
Ubar = 30           # Mean wind going from west to east wind (m/s). (zero by default)
Hbar=5500           # Mean Height (m) for 500 mb surface.

# Part 1. Set constants and domain.
daylen=DAYLEN             #  Total time of integration (in days).
tlen = daylen*24    #  Change to seconds.
Delta_t = DELTA_t*60*60   #  Change to seconds.
nt = tlen/Delta_t         #  Number of time-steps.
t = linspace(0,nt,nt+1)*Delta_t        #  time vector.
time = t/(24*60*60)       #  time in days (for plots).
nspd = (24*60*60)/Delta_t #  time steps per day.
numberoftimes=int(tlen/DELTA_t)  #number of time steps

print(Delta_t)
print(DELTA_t)

#Set grid
nx = NX;ny = NY; nxny = nx*ny                   # Number of points in each direction
print('Grid size, nx=',nx,' ny=',ny);print('Timesteps per day',nspd)

# Calculate the Coriolis parameter and beta parameter
Rearth = 6.371*10**6 #Radius of the Earth (meters).
Omega = 2*pi/(24*60*60) # Angular velocity of the Earth
phi0=45*(pi/180)        # Latitude the calculation is centered on
beta0 = 2*Omega*cos(phi0)/Rearth #Beta parameter
fcor0 = 2*Omega*sin(phi0) #Coiolis parameter

# Calculate the Rossby Radius of Deformation.
grav = 9.81 #gravitational acceleration in ms^-2
L_R = sqrt(grav*Hbar)/fcor0  #Rossby Radius
F = 1/(L_R**2)

#Beta Plane Approximation
fcor0 = 2*Omega*sin(phi0) #Coiolis parameter
y=linspace(111000*-15,111000*15,51)
dfdy45=cos(45*pi/180)*2*Omega/Rearth
fbeta=2*Omega*sin(45*pi/180)+y*dfdy45
fbeta=fbeta[::-1]

fcor0=zeros((nx+1,ny+1))

for i in range(0,51):
    fcor0[i,:]=fbeta[i]



# Specify the domain size and length scale
xlen = 3.33*10**6 # East-West Length of the Domain.
ylen = 3.33*10**6           # North-South Length of the Domain.
Delta_x = xlen/(nx)        # Horizontal grid length
Delta_y = ylen/(ny)        #  Vertical grid length
D_ratio = Delta_y/Delta_x #  Grid length ratio

# Define the grid to have redundant rows east and north.
x = linspace(0,nx,nx+1)*Delta_x
y = linspace(0,ny,ny+1)*Delta_y
(XMESH, YMESH) = meshgrid(x,y);
XX = transpose(XMESH); YY= transpose(YMESH)



# Section 2. Define the Initial Fields.
# The dependent variable is w, the streamfunction. w_0 is the Initial condition.
# w does NOT include the part due to the mean zonal flow.  w is periodic in x and
#symmetric in y. To the constant streamfunction/geopotential field we generate
#a random perturbation in the form of 2-D wave packets

#set seed for same initial conditions every time
seed(a=2)

#randomize seed for different initial conditions
#seed()


Z_0 = zeros((nx+1,ny+1));

# Specify the number of waves in x and y

Nwavex = 1
Nwavey = 1
Nwaves = (2*Nwavex+1)*(2*Nwavey+1)


for kwave in range(-Nwavex,Nwavex+1):
    for lwave in range(-Nwavey,Nwavey+1):
         Amplitude = (2000.0*2*(random()-0.5)) / Nwaves
         phase= 2*pi*(random()-0.5)
         term = cos(2*pi*(kwave*(XX/xlen)+lwave*(YY/ylen)+phase))
         Z_0 = Z_0 + Amplitude*term

#  Add a constant to give typical 500mb values.
Zplus_0 = Z_0 + Hbar

#Add in the zonal mean flow.
Ztotal_0 = Zplus_0 - (fcor0/grav)*Ubar*YY

#Scale the axis and plot the perturbation
XM, YM = XX/(10**3), YY/(10**3)

plt.figure()
CS = plt.contourf(XM, YM, Z_0)
plt.xlabel('x (kilometers)',fontsize=24)
plt.ylabel('y (kilometers)',fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=24)
cb=plt.colorbar()
cb.set_label('Perturbation Height (km)',size=24)
cb.ax.tick_params(labelsize=24)
plt.title('500 mbar Geopotential Perturbation',fontsize=24)

# Plot the field including the mean flow
vecwmin = Ztotal_0.min()
vecwmax = Ztotal_0.max()
vecwmean = (vecwmax+vecwmin)/2
vecwspan = (vecwmax-vecwmin)
vecw = linspace(vecwmean-vecwspan,vecwmean+vecwspan,21)

plt.figure()
plt.contourf(XM, YM, Ztotal_0,vecw)
plt.xlabel('x (kilometers)',fontsize=24)
plt.ylabel('y (kilometers)',fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=24)
cb=plt.colorbar()
cb.set_label('Total Geopotential Height (km)',size=24)
cb.ax.tick_params(labelsize=24)
plt.title('500 mbar Geopotential Height',fontsize=24)

#generate initial streamfunction
w_0=Z_0  #w_0 is perturbation height

#  Add the mean zonal flow
wtotal_0 = w_0 + Hbar - (fcor0/grav)*Ubar*YY;

[XXin,YYin] = meshgrid(x[1:nx],y[1:ny]);
XXin=transpose(XXin)
YYin=transpose(YYin)

R=w_0[0:nx,0:ny]
W_hat = fft2(R)
W_hat_0 = W_hat

#plot the streamfunction as a surface

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(XM, YM, w_0, cmap=cm.coolwarm,
                       linewidth=1)

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.title('Initial Stream Function');
plt.show(block=False)
plt.pause(.000001)


# Section 3. Integrate the BPV Equation in time
#Integrate the BVE in time by leapfrog method. Leapfrog is used to preserve
# energy. Define Q = (Del^2 - F)w. The time derivative is:
# (Q(n+1)-Q(n-1))/(2*Delta_t) and the remaining terms by centered differences:
# R(n) = - ( J(w,Del^2(w)) + beta*(d/dx)w + Ubar*(d/dx)Del^2(w)) at each time
#step. The value of Q at the new time (n+1)*Delta_t is:
#Q(n+1) =  Q(n-1) + 2*Delta_t * R(n). When we have Q(n+1), we have to solve a
#Helmholtz equation to get w(n+1). Then the cycle is repeated.

# Define working arrays to have correct size.
dwdx = zeros((nx+1,ny+1))
dwdy = zeros((nx+1,ny+1))
gradsq = zeros((nx+1,ny+1))
d2wdx2 = zeros((nx+1,ny+1))
d2wdy2 = zeros((nx+1,ny+1))
laplac = zeros((nx+1,ny+1))
dlapdx = zeros((nx+1,ny+1))
dlapdy = zeros((nx+1,ny+1))
wdldx = zeros((nx+1,ny+1))
wdldy = zeros((nx+1,ny+1))
dwdxl = zeros((nx+1,ny+1))
dwdyl = zeros((nx+1,ny+1))
dwdldydx = zeros((nx+1,ny+1))
dwdldxdy = zeros((nx+1,ny+1))
ddwdxldy = zeros((nx+1,ny+1))
ddwdyldx = zeros((nx+1,ny+1))
Jac1 = zeros((nx+1,ny+1))
Jac2 = zeros((nx+1,ny+1))
Jac3 = zeros((nx+1,ny+1))
Jarakawa = zeros((nx+1,ny+1))

# Start of main time-stepping loop
w = w_0;
Energy=[]
Enstrophy=[]
NLsize=[]
CFL_nonlin=[]
R=zeros((nx,ny))
wcenter=[]
zonalwind=[]


plt.figure()
for n in range(1,numberoftimes):

#Take derivatives using finite-difference method, set periodic boundary
#conditions in x and y

#x derivative

    dwdx[1:nx,0:ny+2]= (w[2:nx+1,0:ny+1]-w[0:nx-1,0:ny+1])/(2*Delta_x);
    dwdx[0,0:ny+1] = (w[1,0:ny+1]-w[nx-1,0:ny+1])/(2*Delta_x);
    dwdx[nx,0:ny+1] = dwdx[0,0:ny+1];

# y-derivative of w
    dwdy[0:nx+1,1:ny] = (w[0:nx+1,2:ny+1]-w[0:nx+1,0:ny-1])/(2*Delta_y);

    dwdy[0:nx+1,0] = (w[0:nx+1,1]-w[0:nx+1,ny-1])/(2*Delta_y);
    dwdy[0:nx+1,ny] = dwdy[0:nx+1,0];


# Square of the gradient of w
    gradsq = dwdx**2+dwdy**2;

# Second x-derivative of w
    d2wdx2[1:nx,0:ny+1] = (w[2:nx+1,0:ny+1]+w[0:nx-1,0:ny+1]-2*w[1:nx,0:ny+1])/(Delta_x**2)
    d2wdx2[0,0:ny+1] = (w[1,0:ny+1]+w[nx-1,0:ny+1]-2*w[0,0:ny+1])/(Delta_x**2)
    d2wdx2[nx,0:ny+1] = d2wdx2[0,0:ny+1]

# Second y-derivative of w
    d2wdy2[0:nx+1,1:ny] = (w[0:nx+1,2:ny+1]+w[0:nx+1,0:ny-1]-2*w[0:nx+1,1:ny])/(Delta_y**2)
    d2wdy2[0:nx+1,0] = (w[0:nx+1,1]+w[0:nx+1,ny-1]-2*w[0:nx+1,0])/(Delta_y**2)
    d2wdy2[0:nx+1,ny] = d2wdy2[0:nx+1,0]

    laplac = d2wdx2+d2wdy2;

# x-derivative of laplacian
    dlapdx[1:nx,0:ny+1] = (laplac[2:nx+1,0:ny+1]-laplac[0:nx-1,0:ny+1])/(2*Delta_x)
    dlapdx[0,0:ny+1] = (laplac[1,0:ny+1]-laplac[nx-1,0:ny+1])/(2*Delta_x)
    dlapdx[nx,0:ny+1] = dlapdx[0,0:ny+1]


#y-derivative of laplacian
    dlapdy[0:nx+1,1:ny] = (laplac[0:nx+1,2:ny+1]-laplac[0:nx+1,0:ny-1])/(2*Delta_y)
    dlapdy[0:nx+1,0] = (laplac[0:nx+1,1]-laplac[0:nx+1,ny-1])/(2*Delta_y)
    dlapdy[0:nx+1,ny] = dlapdy[0:nx+1,0]

    Jacobi = dwdx*dlapdy - dwdy*dlapdx

#Compute the Arakawa Jacobian.

    Jac1 = Jacobi;

    wdldx = w*dlapdx
    wdldy = w*dlapdy

    dwdldydx[1:nx,0:ny+1] = (wdldy[2:nx+1,0:ny+1]-wdldy[0:nx-1,0:ny+1])/(2*Delta_x);
    dwdldydx[0,0:ny+1] = (wdldy[1,0:ny+1]-wdldy[nx-1,0:ny+1])/(2*Delta_x);
    dwdldydx[nx,0:ny+1] = dwdldydx[0,0:ny+1];


    dwdldxdy[0:nx+1,1:ny] = (wdldx[0:nx+1,2:ny+1]-wdldx[0:nx+1,0:ny-1])/(2*Delta_y)
    dwdldxdy[0:nx+1,0] = (wdldx[0:nx+1,1]-wdldx[0:nx+1,ny-1])/(2*Delta_y)
    dwdldxdy[0:nx+1,ny] = dwdldxdy[0:nx+1,0]

    Jac2 = dwdldydx - dwdldxdy

    dwdxl = dwdx*laplac
    dwdyl = dwdy*laplac

    ddwdxldy[0:nx+1,1:ny] = (dwdxl[0:nx+1,2:ny+1]-dwdxl[0:nx+1,0:ny-1])/(2*Delta_y)
    ddwdxldy[0:nx+1,0] = (dwdxl[0:nx+1,1]-dwdxl[0:nx+1,ny-1])/(2*Delta_y)
    ddwdxldy[0:nx+1,ny] = ddwdxldy[0:nx+1,0]

    ddwdyldx[1:nx,0:ny+1] = (dwdyl[2:nx+1,0:ny+1]-dwdyl[0:nx-1,0:ny+1])/(2*Delta_x);
    ddwdyldx[0,0:ny+1] = (dwdyl[1,0:ny+1]-dwdyl[nx-1,0:ny+1])/(2*Delta_x)
    ddwdyldx[nx,0:ny+1] = ddwdyldx[0,0:ny+1]

    Jac3 = ddwdxldy - ddwdyldx

    Jarakawa = (1/3)*(Jac1+Jac2+Jac3)

#Compute the function to be stepped forward. The -F*w turn is to dampen out
#small but quickly divergent modes
    Q_n = laplac - F*w

#First time through the loop:
    if n==1:
        Dt = Delta_t/2;
        Q_nm1 = Q_n;

        rmeshvec=linspace(0,nx-1,nx)
        smeshvec=linspace(0,ny-1,ny)
        [rmesh,smesh] = meshgrid(rmeshvec,smeshvec)

        rr = transpose(rmesh)
        ss =transpose(smesh)
        C_rs = 2*(cos(2*pi*rr/nx)-1)/Delta_x**2+2*(cos(2*pi*ss/ny)-1)/Delta_y**2-F

#Calculate the kinetic energy and enstrophy integrals, this allows us to see
#how well energy is being conserved.
    Rgsq=zeros((nx,ny))
    Rgsq[0:nx,0:ny] = gradsq[0:nx,0:ny]
    Rwsq=zeros((nx,ny))
    Rwsq[0:nx,0:ny] = w[0:nx,0:ny]**2
    Energy.append(0.5 * mean(mean(Rgsq+F*Rwsq)))
    Rgsq[0:nx,0:ny] = laplac[0:nx,0:ny]
    Rwsq[0:nx,0:ny] = w[0:nx,0:ny]
    Enstrophy.append(0.5 * mean(mean((Rgsq-F*Rwsq)**2)))
    umax = abs(dwdy).max()
    vmax = abs(dwdx).max()
    maxx=array((umax,vmax))
    VMAX = maxx.max()
    Q_np1 = Q_nm1 - (2*Dt)*((grav/fcor0)*Jacobi + beta0*dwdx + Ubar*dlapdx)

#Section 3.3: Solve the Helmholtz Equation (Del^2-F)w = R.Compute the fft of the
#right hand side (strip off additional row and column).
    R[0:nx,0:ny] = Q_np1[0:nx,0:ny]
    R_hat = fft2(R)

#Compute the transform of the solution
    W_hat = R_hat/C_rs
#Compute the inverse transform to get the solution at (n+1)*Delta_t.
    w_new = real(ifft2(W_hat)) # We assume w is real
    w[0:nx,0:ny] = w_new
    w[nx,0:ny] = w[0,0:ny]     # Fill in additional column at east.
    w[0:nx+1,ny]=w[0:nx+1,0]   # Fill in additional row at north.

    #Add the term for the zonal mean flow.
    wtotal = w + Hbar - (fcor0/grav)*Ubar*YY
    zonalwind.append(-grav*dwdy/fcor0)

#Save an east-west mid cross-section each time-step
#w_section[0:nx+1,n-1] = w[0:nx+1,int(fix(ny/2))]

#Shuffle the fields at the end of each time-step
    Dt = Delta_t
    Q_nm1 = Q_n
    timeestep=n

    #Plot the new 500 mb geopotential height contour
    plt.clf()
    plt.contourf(XM, YM, wtotal,vecw)
    plt.xlabel('x (kilometers)',fontsize=24)
    plt.ylabel('y (kilometers)',fontsize=24)
    plt.tick_params(axis='both', which='major', labelsize=24)
    cb=plt.colorbar()
    cb.set_label('Perturbation Height (km)',size=24)
    cb.ax.tick_params(labelsize=24)
    plt.title('500 mb Geopotential Height and Zonal Wind Contours',fontsize=24)

    #Plot east-west wind strength contours
    plt.draw()
    CS = plt.contour(XM, YM, -grav*dwdy/fcor0,colors='k')
    plt.clabel(CS, inline=1, fontsize=24)
    plt.pause(.00001)

#Plot final height of 500 Millibar pressure level
plt.figure()
plt.contourf(XM, YM, wtotal,vecw)
plt.xlabel('x (kilometers)',fontsize=24)
plt.ylabel('y (kilometers)',fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=24)
cb=plt.colorbar()
cb.set_label('Perturbation Height (km)',size=24)
cb.ax.tick_params(labelsize=24)
plt.title('500 mb Geopotential Height and Zonal Wind Contours after day %s'%(daylen),fontsize=24)


#Plot east-west wind strength contours
plt.draw()
CS = plt.contour(XM, YM, -grav*dwdy/fcor0,colors='k')
plt.clabel(CS, inline=1, fontsize=24)


#Plot the final geopotential anomoly
plt.figure()
CS = plt.contourf(XM, YM, w_0)
plt.xlabel('x (kilometers)',fontsize=24)
plt.ylabel('y (kilometers)',fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=24)
cb=plt.colorbar()
cb.set_label('Perturbation Height (km)',size=24)
cb.ax.tick_params(labelsize=24)
plt.title('500 mbar Geopotential Perturbation after day %s'%(daylen),fontsize=24)

#Plot the final vorticity
plt.figure()
plt.contourf(XM, YM, laplac)
plt.xlabel('x (kilometers)',fontsize=24)
plt.ylabel('y (kilometers)',fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=24)
cb=plt.colorbar()
cb.set_label('Vorticity Anomoly (1/s)',size=24)
cb.ax.tick_params(labelsize=24)
plt.title('500 mb Vorticity Anomoly after day %s'%(daylen),fontsize=24)

#Plot the average east-west wind strength contour over all time-steps
plt.figure()
zonalwindacrosstime = asarray(zonalwind,dtype=float32)
zonalwindavg=mean(zonalwindacrosstime,axis=0)
ZW=plt.contour(XM,YM,zonalwindavg+Ubar)
plt.clabel(ZW, inline=1, fontsize=10)
plt.title('500 mb zonal wind average after day %s'%(daylen),fontsize=24)
plt.show()
