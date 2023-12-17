import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("poisson.txt",skiprows=1)

# 3d graph

from mpl_toolkits.mplot3d import axes3d

ax = plt.figure().add_subplot(projection='3d')

x = data[:,0]
y = data[:,1]
z = data[:,2]
zlim = (1.5, -1.5)
ax.set_zlim(zlim)
ax.plot_trisurf(x, y, z,edgecolor='royalblue', lw=0.5, alpha=0.7, cmap='viridis')
plt.show()