import numpy as np
import matplotlib.pyplot as plt

# Parameter (Anda dapat mengatur nilai-nilai ini sesuai kebutuhan)
m = 1.0  # Massa bola (kg)
g = 9.81  # Percepatan gravitasi (m/s^2)
theta_deg = 30  # Sudut kemiringan (derajat)
v0 = 0.0  # Kecepatan awal (m/s)

# Konversi sudut dari derajat ke radian
theta = np.radians(theta_deg)

# Kondisi awal
x0 = 0.0  # Posisi awal bola di sepanjang sumbu x (m)
y0 = 0.0  # Posisi awal bola di sepanjang sumbu y (m)
vx0 = v0 * np.cos(theta)  # Kecepatan awal dalam arah x (m/s)
vy0 = v0 * np.sin(theta)  # Kecepatan awal dalam arah y (m/s)

# Langkah waktu (delta t)
dt = 0.01  # Delta t (s)

# Jumlah langkah waktu
num_steps = 1000

# Inisialisasi array untuk menyimpan data
t_values = np.zeros(num_steps + 1)
x_values = np.zeros(num_steps + 1)
y_values = np.zeros(num_steps + 1)

# Kondisi awal
t_values[0] = 0.0
x_values[0] = x0
y_values[0] = y0

vyarr = []
vxarr = []
# Metode Euler untuk menghitung posisi bola
for i in range(num_steps):
    t = t_values[i]
    x = x_values[i]
    y = y_values[i]
    
    # Persamaan kecepatan vertikal
    vy = vy0 - g * np.sin(theta) * dt
    
    # persamaan kecepatan horizontal
    vx = vx0 + g * np.cos(theta) * dt

    # Persamaan posisi vertikal
    y_new = y + vy * dt
    
    # Persamaan posisi horizontal
    x_new = x + vx * dt
    
    # Simpan nilai pada array
    t_values[i + 1] = t + dt
    x_values[i + 1] = x_new
    y_values[i + 1] = y_new
    
    # Update kecepatan awal vertikal untuk iterasi berikutnya
    vy0 = vy
    vx0 = vx

z_values = np.zeros(num_steps + 1)

print
# Plot hasil 3d
from matplotlib import animation
# Plot hasil dalam 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ball, =ax.plot(x_values, y_values, z_values, 'o', color='red', markersize=10, zorder = 10)
ax.plot(x_values, y_values, z_values, color='blue', alpha=0.5)
ax.set_xlabel('Posisi Horizontal (m)')
ax.set_ylabel('Posisi Vertikal (m)')
ax.set_zlabel('Posisi Z (m)')
ax.set_title(f'Gerak Bola pada Kemiringan ({theta_deg} derajat, Massa = {m} kg)')

def update(i):
    ball.set_data(x_values[:i], y_values[:i])
    ball.set_3d_properties(0)
    return ball,

ani = animation.FuncAnimation(fig, update, frames=num_steps + 1, interval=10)
# ani.save('test.gif', writer='imagemagick', fps=30)
plt.show()