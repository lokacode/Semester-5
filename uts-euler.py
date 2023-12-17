import numpy as np
import matplotlib.pyplot as plt

# Parameter
m = 1.0
g = 10.0
k = 1.0
v0 = 0.0  # Kecepatan awal
y0 = 0.0  # Posisi awal
t0 = 0.0  # Waktu awal
t_final = 10.0  # Waktu akhir
dt = 0.01  # Langkah waktu

# Inisialisasi array
n_steps = int((t_final - t0) / dt) + 1
t = np.zeros(n_steps)
y = np.zeros(n_steps)
v = np.zeros(n_steps)

# Metode Euler
for i in range(n_steps):
    t[i] = t0 + i * dt
    y[i] = y0
    v[i] = v0

    # Update posisi dan kecepatan
    y0 = y0 + v0 * dt
    v0 = v0 + (g - (k / m) * v0) * dt

print(v[-1])

# invert y
y = -y

# Plot hasil metode Euler
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(t, y)
plt.title('Grafik Posisi (y) terhadap Waktu (Euler)')
plt.xlabel('Waktu (s)')
plt.ylabel('Posisi (m)')

plt.subplot(2, 1, 2)
plt.plot(t, v)
plt.title('Grafik Kecepatan (v) terhadap Waktu (Euler)')
plt.xlabel('Waktu (s)')
plt.ylabel('Kecepatan (m/s)')

plt.tight_layout()
plt.savefig('uts-euler.pdf', dpi=300)
plt.show()
