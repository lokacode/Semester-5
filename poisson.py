def PoissonSolver(u, x, y, nx, ny, eps, Func):

    itmax = 10000 # max no. of iterations

    f = [[0]*(ny+1) for i in range(nx+1)]

    hx = (x[nx]-x[1])/(nx-1); kx = 1e0/(hx*hx) # mesh spacings
    hy = (y[ny]-y[1])/(ny-1); ky = 1e0/(hy*hy)
    kxy = 2e0*(kx + ky)

    for j in range(2,ny): # RHS of PDE for interior points
        for i in range(2,nx): f[i][j] = Func(x[i],y[j])

    for it in range(1,itmax+1): # Gauss-Seidel iteration loop
        err = 0e0
        for j in range(2,ny):
            for i in range(2,nx): # interior mesh points
                uij = (kx*(u[i-1][j] + u[i+1][j]) +
                    ky*(u[i][j-1] + u[i][j+1]) - f[i][j]) / kxy
                eij = 1e0 - uij/u[i][j] if u[i][j] else uij - u[i][j] # local error
                if (abs(eij) > err): err = abs(eij) # maximum error
                u[i][j] = uij

        if (err <= eps): return 0 # convergence test
   
    if (it >= itmax):
        print("PoissonSolver: max. number of iterations exceeded !"); return 0
    return 0

def PoissonXY(u, x, y, nx, ny, eps, Func, CondX, CondY):
    """
    Solves the 2D Poisson equation in Cartesian coordinates on a regular grid
    with (nx x ny) nodes (x[],y[]) using the Gauss-Seidel method. The solution
    u[][] is converged with relative precision eps.
    An error index is returned: 0 - normal execution.
    Calls: Func(x,y) - RHS of Poisson equation; boundary conditions:
    CondX(y,alf_min,bet_min,gam_min,alf_max,bet_max,gam_max)
    CondY(x,alf_min,bet_min,gam_min,alf_max,bet_max,gam_max)
    """
    itmax = 10000 # max no. of iterations
    f = [[0]*(ny+1) for i in range(nx+1)]
    betXmin = [0]*(ny+1); betXmax = [0]*(ny+1)
    gamXmin = [0]*(ny+1); gamXmax = [0]*(ny+1)
    betYmin = [0]*(nx+1); betYmax = [0]*(nx+1)
    gamYmin = [0]*(nx+1); gamYmax = [0]*(nx+1)
    hx = (x[nx]-x[1])/(nx-1); kx = 1e0/(hx*hx) # mesh spacings
    hy = (y[ny]-y[1])/(ny-1); ky = 1e0/(hy*hy)
    kxy = 2e0*(kx + ky)
    
    for j in range(2,ny): # RHS of PDE for interior points
        for i in range(2,nx): f[i][j] = Func(x[i],y[j])
    
    # boundary conditions
    for i in range(1,nx+1): # lower and upper boundaries
        (alf_min,bet_min,gam_min,alf_max,bet_max,gam_max) = CondY(x[i])
        betYmin[i] = bet_min/(alf_min*hy + bet_min)
        gamYmin[i] = gam_min/(alf_min + bet_min/hy)
        betYmax[i] = bet_max/(alf_max*hy + bet_max)
        gamYmax[i] = gam_max/(alf_max + bet_max/hy)
    
    for j in range(2,ny): # left and right boundaries
        (alf_min,bet_min,gam_min,alf_max,bet_max,gam_max) = CondX(y[j])
        betXmin[j] = bet_min/(alf_min*hx + bet_min)
        gamXmin[j] = gam_min/(alf_min + bet_min/hx)
        betXmax[j] = bet_max/(alf_max*hx + bet_max)
        gamXmax[j] = gam_max/(alf_max + bet_max/hx)
    
    for it in range(1,itmax+1): # Gauss-Seidel iteration loop
        err = 0e0
        j=1 # lower boundary
        for i in range(1,nx+1):
            uij = betYmin[i]*u[i][2] + gamYmin[i]
            eij = 1e0 - u[i][j]/uij if uij else uij - u[i][j]
            if (abs(eij) > err): err = abs(eij)
            u[i][j] = uij
        
        for j in range(2,ny):
            i=1 # left boundary
            uij = betXmin[j]*u[i+1][j] + gamXmin[j]
            eij = 1e0 - u[i][j]/uij if uij else uij - u[i][j]
            if (abs(eij) > err): err = abs(eij)
            u[i][j] = uij
        
            for i in range(2,nx): # interior mesh points
                uij = (kx*(u[i-1][j] + u[i+1][j]) +
                    ky*(u[i][j-1] + u[i][j+1]) - f[i][j]) / kxy
                eij = 1e0 - u[i][j]/uij if uij else uij - u[i][j] # local error
                if (abs(eij) > err): err = abs(eij) # maximum error
                u[i][j] = uij
        
            i = nx # right boundary
            uij = betXmax[j]*u[i-1][j] + gamXmax[j]
            eij = 1e0 - u[i][j]/uij if uij else uij - u[i][j]
            if (abs(eij) > err): err = abs(eij)
            u[i][j] = uij
        
        j = ny # upper boundary
        for i in range(1,nx+1):
            uij = betYmax[i]*u[i][ny-1] + gamYmax[i]
            eij = 1e0 - u[i][j]/uij if uij else uij - u[i][j]
            if (abs(eij) > err): err = abs(eij)
            u[i][j] = uij
        
        if (err <= eps): break # convergence test
    if (it >= itmax):
        print("PoissonXY: max. number of iterations exceeded !"); return 1
    return 0

from math import cos, pi


def Func(x, y): 
    return cos(x + y) - cos(x - y)

def CondX(y): 
    alf_min = 1e0
    bet_min = 0e0
    gam_min = 0e0
    alf_max = 1e0
    bet_max = 0e0
    gam_max = 0e0
    return (alf_min, bet_min, gam_min, alf_max, bet_max, gam_max)

def CondY(x): 
    alf_min = 1e0
    bet_min = 0e0
    gam_min = 0e0
    alf_max = 1e0
    bet_max = 0e0
    gam_max = 0e0
    return (alf_min, bet_min, gam_min, alf_max, bet_max, gam_max)

# main
xmin = -pi
xmax = pi
ymin = -pi
ymax = pi
nx = 51
ny = 51
eps = 1e-5

u = [[0]*(ny+1) for i in range(nx+1)]
x = [0]*(nx+1)
y = [0]*(ny+1)

hx = (xmax-xmin)/(nx-1)
for i in range(1, nx+1):
    x[i] = xmin + (i-1)*hx 

hy = (ymax-ymin)/(ny-1)
for j in range(1, ny+1):
    y[j] = ymin + (j-1)*hy 

for j in range(1, ny+1):
    for i in range(1, nx+1):
        u[i][j] = 0e0

PoissonXY(u, x, y, nx, ny, eps, Func, CondX, CondY)

out = open("Poisson.txt", "w")
out.write(" x y u\n")

for j in range(1, ny+1):
    for i in range(1, nx+1):
        out.write(("{0:10.5f}{1:10.5f}{2:14.5e}\n").format(x[i], y[j], u[i][j]))

out.close()

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("poisson.txt",skiprows=1)

# 3d graph
ax = plt.figure().add_subplot(projection='3d')

x = data[:,0]
y = data[:,1]
z = data[:,2]
zlim = (1.5, -1.5)
ax.set_zlim(zlim)
ax.plot_trisurf(x, y, z,edgecolor='royalblue', lw=0.5, alpha=0.7, cmap='viridis')
plt.savefig('poissond.png', dpi=300, bbox_inches='tight')
plt.show()