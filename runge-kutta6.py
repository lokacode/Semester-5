def RungeKutta(t, ht, y, n, Func):
    """
    Propagates the solution y of a system of 1st order ODEs
    y'[i] = f[i](t, y[]), i = 1..n
    from t to t+ht using the 4th order Runge-Kutta method.
    Calls: Func(t, y, f) - RHS of ODEs
    """
    f1 = [0] * (n + 1)
    f2 = [0] * (n + 1)
    f3 = [0] * (n + 1)
    f4 = [0] * (n + 1)
    f5 = [0] * (n + 1)
    f6 = [0] * (n + 1)
    yt = [0] * (n + 1)  # predicted solution
    ht2 = ht / 2e0

    # Calculate the RHS at t
    Func(t, y, f1)

    for i in range(1, n + 1):
        yt[i] = y[i] + ht/4 * f1[i]
    Func(t + ht/4, yt, f2)

    for i in range(1, n + 1):
        yt[i] = y[i] + (ht/8 * f1[i]) + (ht/8 * f2[i])
    Func(t + ht/4, yt, f3)

    for i in range(1, n + 1):
        yt[i] = y[i] - ht/2 * f2[i] + ht * f3[i]
    Func(t + ht/2, yt, f4)

    for i in range(1, n + 1):
        yt[i] = y[i] + (3 * ht/16 * f1[i]) + (9 * ht/16 * f4[i])
    Func(t + ht*3/4, yt, f5)

    for i in range(1, n + 1):
        yt[i] = y[i] - (3 * ht/7 * f1[i]) + (2 * ht/7 * f2[i]) + (12 * ht/7 * f3[i]) - (12 * ht/7 * f4[i]) + (8 * ht/7 * f5[i])
    Func(t + ht, yt, f6)

    # Propagate the solution
    for i in range(1, n + 1):
        y[i] += ht/90 * (7 * f1[i] + 32 * f3[i] + 12 * f4[i] + 32 * f5[i] + 7 * f6[i])
# Angular motion of a nonlinear pendulum by the Runge-Kutta method
# u" = -g/l * sin(u) - k * u’, u(0) = u0, u’(0) = u0’
from math import *
import numpy as np

g = 9.81e0 # gravitational acceleration
def Func(t, u, f): # RHS of 1st order ODEs
    f[1] = u[2] # u[1] = u, u[2] = u’
    f[2] = -g/l * sin(u[1]) - k * u[2]

# main
l = 1e0 # pendulum length
k = 0e0 # velocity coefficient
# initial displacement
u0 = [0.1, pi/2, 3.1] # initial displacement
du0 = 0e0 # initial derivative
tmax = 20e0 # time span
ht = 0.001e0 # time step size

for i in range(len(u0)):
    n=2 # number of 1st order ODEs
    u = [0]*(n+1) # solution components

    out = open(f"pendulum{i}.txt","w") # open output file
    out.write(" t u du\n")

    t = 0e0
    u[1] = u0[i]; u[2] = du0 # initial values
    out.write(("{0:10.5f}{1:10.5f}{2:10.5f}\n").format(t,u[1],u[2]))

    nT = 0 # number of half-periods
    t1 = t2 = 0e0 # bounding solution zeros
    us = u[1] # save solution

    while (t+ht <= tmax): # propagation loop
        RungeKutta(t,ht,u,n,Func)
        t += ht
        if (u[1]*us < 0e0): # count solution passages through zero
            if (t1 == 0): t1 = t # initial zero
            else: t2 = t; nT += 1 # final zero
        us = u[1] # save solution
        out.write(("{0:10.5f}{1:10.5f}{2:10.5f}\n").format(t,u[1],u[2]))

    T = 2e0*(t2-t1) / nT # calculated period
    T0 = 2e0*pi*sqrt(l/g) # harmonic period
    print("u0 = {0:7.5f} T/T0 = {1:7.5f}".format(u0[i],T/T0))
    out.close()

u0 = [0.1, 0.6, 1.1, 1.6, 2.1, 2.52, 2.8, 3, 3.1]
pe = []

# Runge-Kutta
for i in range (len(u0)):
    t = 0e0
    u[1] = u0[i]; u[2] = du0 # initial values
    nT = 0 # number of half-periods
    t1 = t2 = 0e0 # bounding solution zeros
    us = u[1] # save solution
    while (t+ht <= tmax): # propagation loop
        RungeKutta(t,ht,u,n,Func)
        t += ht
        if (u[1]*us < 0e0): # count solution passages through zero
            if (t1 == 0): t1 = t # initial zero
            else: t2 = t; nT += 1 # final zero
        us = u[1] # save solution

    T = 2e0*(t2-t1) / nT # calculated period
    T0 = 2e0*pi*sqrt(l/g) # harmonic period
    pe.append(T/T0)

u1 = np.linspace(0.1,3.1,100)
pd = []

# theory
for i in range (len(u1)):
    t = 0e0
    u[1] = u1[i]; u[2] = du0 # initial values
    nT = 0 # number of half-periods
    t1 = t2 = 0e0 # bounding solution zeros
    us = u[1] # save solution
    while (t+ht <= tmax): # propagation loop
        RungeKutta(t,ht,u,n,Func)
        t += ht
        if (u[1]*us < 0e0): # count solution passages through zero
            if (t1 == 0): t1 = t # initial zero
            else: t2 = t; nT += 1 # final zero
        us = u[1] # save solution

    T = 2e0*(t2-t1) / nT # calculated period
    T0 = 2e0*pi*sqrt(l/g) # harmonic period
    pd.append(T/T0)

# Plotting
import matplotlib.pyplot as plt

# plt.xkcd()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

data = np.loadtxt("pendulum0.txt",skiprows=1)
data1= np.loadtxt("pendulum1.txt",skiprows=1)
data2= np.loadtxt("pendulum2.txt",skiprows=1)

ax1.plot(data[:,0],data[:,1],label='u')
ax1.plot(data1[:,0],data1[:,1],label='u1')
ax1.plot(data2[:,0],data2[:,1],label='u2')
ax1.set_xlabel('t (s)')
ax1.set_ylabel('u (rad)')
ax1.set_xlim(0,20)
ax1.set_ylim(-4,4)
ax1.tick_params(direction='in', top=True, right=True)
# ax1.legend()

# scatter plot
ax2.scatter(u0,pe,label='Runge-Kutta', edgecolors='black', facecolors='none')
# line plot
ax2.plot(u1,pd,label='theory', color='darkblue', linewidth=1)
ax2.set_xlabel('u0 (rad)')
ax2.set_ylabel('T/T0')
ax2.set_xlim(0,3.5)
ax2.set_ylim(0.5,3.5)
ax2.tick_params(direction='in', top=True, right=True)
ax2.legend()
plt.tight_layout()
# plt.savefig('pendulum.png', dpi=300)
plt.show()