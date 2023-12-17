import numpy as np
import matplotlib.pyplot as plt

# Parameter
m = 1.0
g = 10.0
k = 1
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

# Metode Runge-Kutta orde 4
for i in range(n_steps):
    t[i] = t0 + i * dt
    y[i] = y0
    v[i] = v0

    # K1, K2, K3, K4
    K1v = (g - (k / m) * v0) * dt
    K1y = v0 * dt

    K2v = (g - (k / m) * (v0 + 0.5 * K1v)) * dt
    K2y = (v0 + 0.5 * K1v) * dt

    K3v = (g - (k / m) * (v0 + 0.5 * K2v)) * dt
    K3y = (v0 + 0.5 * K2v) * dt

    K4v = (g - (k / m) * (v0 + K3v)) * dt
    K4y = (v0 + K3v) * dt

    # Update posisi dan kecepatan
    y0 = y0 + (K1y + 2 * K2y + 2 * K3y + K4y) / 6
    v0 = v0 + (K1v + 2 * K2v + 2 * K3v + K4v) / 6

print(v[-1])

# invert y
y = -y

# export hasil ke excel
# import pandas as pd
# df = pd.DataFrame({'t': t, 'y': y, 'v': v})
# df.to_excel('soal2021-rk4.xlsx', index=False)

# Plot hasil metode Runge-Kutta
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(t, y)
plt.title('Grafik Posisi (y) terhadap Waktu (Runge-Kutta)')
plt.xlabel('Waktu (s)')
plt.ylabel('Posisi (m)')

plt.subplot(2, 1, 2)
plt.plot(t, v)
plt.title('Grafik Kecepatan (v) terhadap Waktu (Runge-Kutta)')
plt.xlabel('Waktu (s)')
plt.ylabel('Kecepatan (m/s)')

plt.tight_layout()
plt.savefig('uts-rk4.pdf', dpi=300)
plt.show()
