# Arden Ahmad
# 121110018

# Reflexion/transmission of quantum wave packet using Gauss-Seidel solver
from math import *
def Pot(x, a, V0): # Potential
    return V0 if fabs(x) <= 0.5e0*a else 0e0
def PropagQGS(Psi, Chi, V, nx, hx, ht):
#----------------------------------------------------------------------------
# Propagates the solution PSI = Psi + i Chi of the 1D Schrodinger equation
# i d/dt PSI(x,t) = [-(1/2) d2/dx2 - V(x)] PSI(x,t),
# over the time interval ht. Uses the Crank-Nicolson scheme on a grid with
# nx nodes and spacing hx and solves the discretized system by Gauss-Seidel
# iterations. Uses real arithmetic.
#----------------------------------------------------------------------------
    eps = 1e-6 # relative tolerance
    itmax = 100 # max no. of iterations
    Psi0 = [0]*(nx+1); Chi0 = [0]*(nx+1); W = [0]*(nx+1)
    lam = ht/(4e0*hx*hx)
    for i in range(2,nx): # explicit 0th-order approximations
        W[i] = 0.5e0*ht*V[i] + 2e0*lam
        Psi0[i] = Psi[i] - lam*(Chi[i-1] + Chi[i+1]) + W[i]*Chi[i]
        Chi0[i] = Chi[i] + lam*(Psi[i-1] + Psi[i+1]) - W[i]*Psi[i]
    for it in range(1,itmax+1): # Gauss-Seidel iteration loop
        err = 0e0
        for i in range (2,nx):
            Psii = Psi0[i] - lam*(Chi[i-1] + Chi[i+1]) + W[i]*Chi[i]
            Chii = Chi0[i] + lam*(Psi[i-1] + Psi[i+1]) - W[i]*Psii
            # local error estimate based on probability density
            erri = fabs((Psii*Psii+Chii*Chii) - (Psi[i]*Psi[i]+Chi[i]*Chi[i]))
            if (erri > err): err = erri # maximum error estimate
            Psi[i] = Psii
            Chi[i] = Chii
        if (err <= eps): break # convergence test
    if (it > itmax):
        print("PropagQGS: max. number of iterations exceeded !"); return 1
    return 0
#============================================================================
def InitR(Psi, Chi, x, nx, x0, sig, k0):
#----------------------------------------------------------------------------
# Initial Gaussian wave packet
# Psi(x,0) = 1/sqrt(sqrt(2*pi)*sig) * exp[-(x-x0)^2/(4*sig^2)] * exp(ikx)
# x0 - position of center, sig - half-width, k0 - average wave number
#----------------------------------------------------------------------------
    a = 1e0/sqrt(sqrt(2e0*pi)*sig)
    b =-1e0/(4e0*sig*sig)
    for i in range(1,nx+1):
        dx = x[i] - x0
        f=a * exp(b*dx*dx)
        if (f < 1e-10): f = 0e0
        Psi[i] = f * cos(k0*x[i]); Chi[i] = f * sin(k0*x[i])
#============================================================================
def ProbDensR(Psi, Chi, Psi2, nx, hx):
#----------------------------------------------------------------------------
# Calculates the probability density Psi2[] of the wave function Psi[]
#----------------------------------------------------------------------------
    for i in range(1,nx+1): # unnormalized probability density
        Psi2[i] = Psi[i]*Psi[i] + Chi[i]*Chi[i]
        if (Psi2[i] <= 1e-10): Psi2[i] = 0e0
    PsiNorm = 0.5e0*(Psi2[1] + Psi2[nx]) # integral by trapezoidal rule
    for i in range(2,nx): PsiNorm += Psi2[i]
    PsiNorm *= hx
    for i in range(1,nx+1): Psi2[i] /= PsiNorm # normalized prob. density
    return PsiNorm
# main
a = 5e0 # width of potential barrier
V0 = 30e0 # height of potential barrier
x0 = -20e0 # initial position of wave packet
sig = 1e0 # half-width of packet
k0 = 10e0 # average wave number of packet
xmax = 100e0 # maximum x
hx = 1e-1 # spatial step size
tmax = 5e0 # maximum propagation time
ht = 1.25e-3 # time step
nout = 40 # output every nout steps
nx = 2*(int)(xmax/hx + 0.5) + 1 # odd number of spatial nodes
nt = (int)(tmax/ht + 0.5) # number of time steps
nx2 = int(nx/2)
Psi = [0]*(nx+1) # real part of wave function
Chi = [0]*(nx+1) # imag part of wave function
Psi2 = [0]*(nx+1) # probability density
V = [0]*(nx+1) # potential
x = [0]*(nx+1) # spatial mesh
for i in range(1,nx+1): # tabulate spatial mesh and potential
    x[i] = (i-nx2-1)*hx
    V[i] = Pot(x[i],a,V0)
InitR(Psi,Chi,x,nx,x0,sig,k0) # initial wave packet
for it in range(1,nt+1): # time loop
    t = it*ht
    PropagQGS(Psi,Chi,V,nx,hx,ht) # propagate by Gauss-Seidel solver
    PsiNorm = ProbDensR(Psi,Chi,Psi2,nx,hx) # probability density
    if (it % nout == 0 or it == nt): # output every nout steps
        fname = "scatter_{0:4.2f}.txt".format(t)
        out = open(fname,"w")
        out.write("t = {0:4.2f}\n".format(t))
        out.write(" x V PsiR PsiI Psi2\n")
        for i in range(1,nx+1):
            out.write("{0:10.5f}{1:10.5f}{2:10.5f}{3:10.5f}{4:10.5f}\n".\
                format(x[i],V[i],Psi[i],Chi[i],Psi2[i]))
        out.close

# animation
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# Create a function to update the plot for each frame
def update(frame):
    plt.clf()  # Clear the previous frame
    fname = "scatter_{0:4.2f}.txt".format((frame + 1) * 0.05)
    data = np.loadtxt(fname, skiprows=2)
    x = data[:, 0]
    Psi2 = data[:, 4]
    psir = data[:, 1]
    plt.plot(x, Psi2, label='Psi2')
    plt.plot(x, psir, label='psir')
    plt.xlabel('x')
    plt.xlim(-30, 30)
    plt.ylim(-0.1, 0.5)
    plt.ylabel('Psi2')
    plt.title('t = {0:4.2f}'.format((frame + 1) * 0.05))
    # plt.legend(loc='upper right')

# Create a figure
fig, ax = plt.subplots()

# Create the animation
ani = FuncAnimation(fig, update, frames=100, interval=100)

# Save the animation as a video (e.g., MP4)
ani.save('scatter_animation.mp4', writer='ffmpeg', dpi=300)

# Display the animation (optional)
plt.show()
