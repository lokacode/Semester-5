#============================================================================
def Euler(t, ht, y, n, Func):
#----------------------------------------------------------------------------
# Propagates the solution y of a system of 1st order ODEs
# y’[i] = f[i](t,y[]), i = 1..n
# from t to t+ht using Euler’s method
# Calls: Func(t, y, f) - RHS of ODEs
#----------------------------------------------------------------------------
    f = [0]*(n+1) # RHS of ODEs
    Func(t,y,f) # get RHS of ODEs
    for i in range(1,n+1): y[i] += ht * f[i] # propagate solution
#============================================================================
def EulerPC(t, ht, y, n, Func):
#----------------------------------------------------------------------------
# Propagates the solution y of a system of 1st order ODEs
# y’[i] = f[i](t,y[]), i = 1..n
# from t to t+ht using Euler’s predictor-corrector method
# Calls: Func(t, y, f) - RHS of ODEs
#----------------------------------------------------------------------------
    f1 = [0]*(n+1); f2 = [0]*(n+1) # RHS of ODEs
    yt = [0]*(n+1) # predicted solution
    Func(t,y,f1) # RHS of ODEs at t
    for i in range(1,n+1): yt[i] = y[i] + ht * f1[i] # predictor
    Func(t+ht,yt,f2) # RHS of ODEs at t+ht
    ht2 = ht/2e0
    for i in range(1,n+1): y[i] += ht2 * (f1[i] + f2[i]) # corrector

# Solves a Cauchy problem for a 2nd order ODE by Euler’s method
# y" + y = 0, y(0) = y0, y’(0) = y0’
# Equivalent problem: y[1] = y, y[2] = y’
# y[1]’ = y[2], y[1](0) = y0
# y[2]’ = -y[1], y[2](0) = y0’
#----------------------------------------------------------------------------
from math import *
def Func(t, y, f): # Right-hand sides of ODEs
    f[1] = y[2]
    f[2] = -y[1]
# main
y0 = 0e0; dy0 = 1e0 # initial values => y(t) = sin(t)
tmax = 100e0 # time span
ht = 0.05e0 # step size
n=2 # number of 1st order ODEs
y = [0]*(n+1) # solution components

# euler
out = open("ode.txt","w") # open output file
out.write(" t y1 y2 check\n")
t = 0e0
y[1] = y0; y[2] = dy0 # initial values
out.write(("{0:10.5f}{1:10.5f}{2:10.5f}{3:10.5f}\n").format(t,y[1],y[2],y[1]*y[1]+y[2]*y[2]))
while (t+ht <= tmax): # propagation loop
    Euler(t,ht,y,n,Func)
    t += ht
    out.write(("{0:10.5f}{1:10.5f}{2:10.5f}{3:10.5f}\n").format(t,y[1],y[2],y[1]*y[1]+y[2]*y[2]))
out.close()

# euler predictor-corrector
out = open("odepc.txt","w") # open output file
out.write(" t y1 y2 check\n")
t = 0e0
y[1] = y0; y[2] = dy0 # initial values
out.write(("{0:10.5f}{1:10.5f}{2:10.5f}{3:10.5f}\n").format(t,y[1],y[2],y[1]*y[1]+y[2]*y[2]))
while (t+ht <= tmax): # propagation loop
    EulerPC(t,ht,y,n,Func)
    t += ht
    out.write(("{0:10.5f}{1:10.5f}{2:10.5f}{3:10.5f}\n").format(t,y[1],y[2],y[1]*y[1]+y[2]*y[2]))
out.close()

import matplotlib.pyplot as plt
import numpy as np
# plt.style.use('seaborn-poster')

title_font = {'fontname': 'Times New Roman', 'fontsize': 18}
label_font = {'fontname': 'Times New Roman', 'fontsize': 14}

data1 = np.loadtxt("ode.txt", skiprows=1)
data2 = np.loadtxt("odepc.txt", skiprows=1)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 12))

# Plot for data1
ax1.plot(data1[:, 1], data1[:, 2], 'r-')
ax1.set_xlabel('y1', **label_font)
ax1.set_ylabel('y2', **label_font)
ax1.set_title('Euler Method', **title_font)
ax1.tick_params(axis='both', which='major', labelsize=12)

# Plot for data2
ax2.plot(data2[:, 1], data2[:, 2], '-')
ax2.set_xlabel('y1', **label_font)
ax2.set_ylabel('y2', **label_font)
ax2.set_title('Euler Predictor-Corrector', **title_font)
ax2.tick_params(axis='both', which='major', labelsize=12)

# Plot for data1 sin(t)
ax3.plot(data1[:, 0], data1[:, 1], 'r-')
ax3.set_xlabel('t', **label_font)
ax3.set_ylabel('y1', **label_font)
ax3.set_title('Euler Method', **title_font)
ax3.tick_params(axis='both', which='major', labelsize=12)
ax3.set_xlim(0, 100)
ax3.axhline(y = 0, color='dimgrey', linestyle='-', linewidth=1)
ax3.grid(which='major', linestyle='--', linewidth='0.5', color='dimgrey')

# Plot for data2 sin(t)
ax4.plot(data2[:, 0], data2[:, 1], '-')
ax4.set_xlabel('t', **label_font)
ax4.set_ylabel('y1', **label_font)
ax4.set_title('Euler Predictor-Corrector', **title_font)
ax4.tick_params(axis='both', which='major', labelsize=12)
ax4.set_ylim(-2, 2)
ax4.set_xlim(0, 100)
ax4.axhline(y = 0, color='dimgrey', linestyle='-', linewidth=1)
ax4.grid(which='major', linestyle='--', linewidth='0.5', color='dimgrey')

plt.tight_layout()
# plt.savefig('ode1.png', dpi=300)
plt.show()