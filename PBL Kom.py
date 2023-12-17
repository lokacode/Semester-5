import numpy as np
import matplotlib.pyplot as plt

# Parameter objek massif
M = 5.972e24  # Massa Bumi dalam kilogram
R = 6371e3  # Jari-jari Bumi dalam meter

# Parameter grid
N = 100  # Jumlah titik dalam satu dimensi
x_max, y_max = 2 * R, 2 * R  # Rentang koordinat x dan y
dx = x_max / N  # Ukuran langkah
dy = y_max / N  # Ukuran langkah

# Membuat grid koordinat
x = np.linspace(-x_max, x_max, N)
y = np.linspace(-y_max, y_max, N)
X, Y = np.meshgrid(x, y)

# Inisialisasi medan gravitasi
phi = np.zeros((N, N))

# Mendefinisikan fungsi yang menghitung medan gravitasi
def calculate_gravitational_field(phi, M, R, dx, dy):
    for i in range(1, N - 1):
        for j in range(1, N - 1):
            r = np.sqrt(X[i, j] ** 2 + Y[i, j] ** 2)
            phi[i, j] = -G * M / r

# Konstanta gravitasi
G = 6.67430e-11  # Konstanta gravitasi universal Newton

# Menghitung medan gravitasi
calculate_gravitational_field(phi, M, R, dx, dy)

# Membuat plot 2D distribusi medan gravitasi
plt.figure(figsize=(10, 8))
plt.contourf(X, Y, phi, levels=50, cmap="viridis")
plt.colorbar(label="Medan Gravitasi (N/kg)")
plt.xlabel("Koordinat X (meter)")
plt.ylabel("Koordinat Y (meter)")
plt.title("Distribusi Medan Gravitasi di Sekitar Objek Massif (Bumi)")
plt.axis('equal')
plt.show()
