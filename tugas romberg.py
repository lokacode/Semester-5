# Arden Ahmad
# 121110018

from math import *
import matplotlib.pyplot as plt
import numpy as np

#============================================================================
def qRomberg(Func, a, b, eps = 1e-6):
#----------------------------------------------------------------------------
# Integrates function Func on interval [a,b] with relative precision eps
# using the adaptive Romberg method
#----------------------------------------------------------------------------
    kmax = 30 # max. no. of step halving iterations
    r1 = [0]*(kmax+1) # two consecutive lines
    r2 = [0]*(kmax+1) # from the method table
    h = b-a; n = 1

    r1[0] = 0.5*h*(Func(a) + Func(b)) # initial approximation
    for k in range(1,kmax+1): # step halving loop
        sumf = 0e0
        for i in range(1,n+1): sumf += Func(a+(i-0.5)*h)
        r2[0] = 0.5*(r1[0] + h*sumf) # trapezoid formula
        f = 1e0
        for j in range(1,k+1): # increase quadrature order
            f *= 4
            r2[j] = (f*r2[j-1] - r1[j-1])/(f-1) # new approximation

        if (k > 1): # convergence check
            if (fabs(r2[k]-r1[k-1]) <= eps*fabs(r2[k])): break
            if (fabs(r2[k]) <= eps and fabs(r2[k]) <= fabs(r2[k]-r1[k-1])):break
        h *= 0.5; n *= 2 # halve integration step
        for j in range(0,k+1): r1[j] = r2[j] # shift table lines

    if (k >= kmax):
        print("qRomberg: max. no. of iterations exceeded !")
        k -= 1

    return r2[k]

# constants
g = 9.81 # gravitational acceleration (m/s^2)
m = 80  # mass (kg)
c = 10 # drag coefficient (kg/s)


# function to be integrated
def Func(t):
    return g*m/c*(1 - exp(-(c/m)*t))

# integration interval
a = 0e0; b = 8e0

# Romberg integration
eps = 1e-6 # relative precision
I = qRomberg(Func, a, b, eps)
print(f"Hasil = {I:.6f}")

# area under the curve
x = np.linspace(-2,10,1000)
x1 = np.linspace(a,b,1000)
y = [Func(t) for t in x]
y1 = [Func(t) for t in x1]

# plot the function
plt.plot(x,y)
plt.fill_between(x1,y1,0,facecolor='cyan',alpha=0.5)
# plt.axis([a,8,0,100])
plt.xlabel('t (s)'); plt.ylabel('v (m/s)')
plt.title(r'$v = \frac{g \cdot m}{c} \left(1 - e^{-(c/m) \cdot t}\right)$', fontsize=14)
plt.grid(axis='x')
# x and y axis line
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
# remove down and left ticks
plt.gca().spines['bottom'].set_position(('data',0))
plt.gca().spines['left'].set_position(('data',0))
# remove up and right ticks
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
# text result
plt.text(2, 50, f"Luas (jarak) = {I:.6f}", fontsize=12)
plt.tight_layout()
plt.savefig('tugas romberg.png', dpi=300)
plt.show()
