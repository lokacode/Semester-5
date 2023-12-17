import numpy as np
import matplotlib.pyplot as plt

# Membuat data
x = np.linspace(-3, 3, 100)
y = np.linspace(-3, 3, 100)
X, Y = np.meshgrid(x, y)
Z = np.sin(X) * np.sin(Y)

# Membuat plot 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')

# Menambahkan label sumbu
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_zlim(1.5, -1.5)
ax.set_ylim(3, -3)
ax.set_xlim(3, -3)
# Menampilkan grafik
plt.show()