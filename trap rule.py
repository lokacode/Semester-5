# Evaluates an integral using the trapezoidal rule
from math import *
from time import time
import numpy as np

def Func(x): return (x*x*x) * exp(-x)
#============================================================================
def qTrapz(Func, a, b, n):
#----------------------------------------------------------------------------
# Integrates function Func on interval [a,b] using the trapezoidal rule
# with n integration points
#----------------------------------------------------------------------------
    h = (b-a)/(n-1)
    print(f"h = {h}")
    s = 0.5*(Func(a) + Func(b))
    for i in range(1,n-1): s += Func(a+i*h)
    return h*s


#============================================================================
def qSimpson(Func, a, b, n):
#----------------------------------------------------------------------------
# Integrates function Func on interval [a,b] using Simpson’s rule with n
# (odd) integration points
#----------------------------------------------------------------------------
    if (n % 2 == 0): n += 1 # increment n if even
    h = (b-a)/(n-1)
    s1 = s2 = 0e0
    for i in range(2,n-2,2): s1 += Func(a + i*h) # odd-index sum
    for i in range(1,n-1,2): s2 += Func(a + i*h) # even-index sum
    return (h/3)*(Func(a) + 4*s2 + 2*s1 + Func(b))

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


analitik = 0.11392894125692

# integration interval and number of integration points
a = 0e0; b = 1e0; n = 100
time1 = time()
print("\nI trapezoidal = ",qTrapz(Func,a,b,n))
time2 = time()
print("Elapsed time: ",time2-time1," sec")

time1 = time()
print("\nI Romberg = ",qRomberg(Func,a,b))
time2 = time()
print("Elapsed time: ",time2-time1," sec")

time1 = time()
print("\nI simpson = ",qSimpson(Func,a,b,n))
time2 = time()
print("Elapsed time: ",time2-time1," sec")


print("\nAnalitik = ",analitik)