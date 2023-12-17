import numpy as np
import matplotlib.pyplot as plt

# Konstanta-konstanta yang diperlukan
G = 6.67430e-11  # Konstanta gravitasi universal (m^3/kg/s^2)
mass_bumi = 5.972e24  # Massa Bumi (kg)
radius_bumi = 6371e3  # Radius Bumi (m)

# Membuat grid
grid = 500  # Resolusi grid (banyak titik)
x = np.linspace(-2*radius_bumi, 2*radius_bumi, grid)
y = np.linspace(-2*radius_bumi, 2*radius_bumi, grid)
X, Y = np.meshgrid(x, y)
phi = np.zeros_like(X)

print(f'Δr = {x[1] - x[0]} m') # Ukuran langkah (Δr)
print(f'r0 = {x[0]} m') # Titik awal (r0)
print(f'rmax = {x[-1]} m') # Titik akhir (rmax)
print(f'N = {grid**2} titik') # Jumlah titik (N)

# Menghitung medan gravitasi (N/kg) dengan metode beda hingga
for i in range(grid):
    for j in range(grid):
        r = np.sqrt(x[i]**2 + y[j]**2)  # Jarak dari pusat Bumi
        if r >= radius_bumi:
            phi[i, j] = G * mass_bumi / r**2  # Medan gravitasi di luar Bumi
        else:
            phi[i, j] = G * mass_bumi / radius_bumi**2  # Medan gravitasi di permukaan Bumi

# Plot Medan Gravitasi
from matplotlib.patches import Arrow
arrow = Arrow(0, 0, radius_bumi, 0, color='red', width= 50, label='Radius Bumi')
plt.figure(figsize=(9, 7))
plt.imshow(phi, extent=(-2*radius_bumi, 2*radius_bumi, -2*radius_bumi, 2*radius_bumi), cmap='viridis', origin='lower')
plt.colorbar(label='Medan Gravitasi (N/kg)')
plt.title('Medan Gravitasi di Sekitar Bumi')
plt.xlabel('Jarak (m)')
plt.ylabel('Jarak (m)')
plt.gca().set_aspect('equal', adjustable='box')
plt.gca().add_patch(arrow)
plt.legend(loc='upper left')
plt.tight_layout()
# plt.savefig('medan_gravitasi.png', dpi=300)
plt.show()