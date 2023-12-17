import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# Konstanta
WIDTH, HEIGHT = 80, 40
CELL_SIZE = 1

# Variabel untuk simulasi gelombang air dan angin
damping = 0.99

# Variasi kecepatan, sudut, dan tekanan angin
wind_speeds = np.linspace(0.05, 0.5, WIDTH)
wind_angles = np.linspace(0, 2*np.pi, HEIGHT)
pressure_map = np.random.uniform(0.1, 0.5, size=(WIDTH, HEIGHT))

# Matriks tinggi air
height_map = np.zeros((WIDTH, HEIGHT), dtype=float)
velocity_map = np.zeros_like(height_map, dtype=float)

# Fungsi untuk mensimulasikan gelombang air dengan pengaruh angin
def simulate_water():
    global height_map, velocity_map
    new_height_map = height_map.copy()
    new_velocity_map = velocity_map.copy()

    for i in range(1, WIDTH - 1):
        for j in range(1, HEIGHT - 1):
            # Perhitungan tinggi air dan kecepatan dengan mempertimbangkan angin
            wind_speed = wind_speeds[i]
            wind_angle = wind_angles[j]
            wind_force = pressure_map[i, j] * wind_speed

            new_height_map[i, j] = (
                height_map[i, j] +
                (height_map[i - 1, j] + height_map[i + 1, j] + height_map[i, j - 1] + height_map[i, j + 1]) / 4 -
                new_height_map[i, j]
            ) * damping - wind_force * np.cos(wind_angle)

            new_velocity_map[i, j] = (
                velocity_map[i, j] +
                (velocity_map[i - 1, j] + velocity_map[i + 1, j] + velocity_map[i, j - 1] + velocity_map[i, j + 1]) / 4 -
                new_velocity_map[i, j]
            ) * damping

    height_map = new_height_map
    velocity_map = new_velocity_map

# Inisialisasi plot 3D
fig = plt.figure(figsize=(8, 8), tight_layout=True)
ax = fig.add_subplot(111, projection='3d')

# Fungsi animasi
def update(frame):
    simulate_water()
    ax.clear()
    x, y = np.meshgrid(wind_speeds, wind_angles)
    surface = ax.plot_surface(x, y, height_map.T, cmap='Blues', rstride=1, cstride=1, alpha=0.8)
    ax.set_title('Simulasi Gelombang Air 3D')
    # set y ticks label dalam derajat
    ax.set_yticks(np.linspace(0, 2*np.pi, 5))
    ax.set_yticklabels([f'{int(np.degrees(angle))}°' for angle in np.linspace(0, 2*np.pi, 5)])
    ax.set_xlabel('Kecepatan Angin')
    ax.set_ylabel('Sudut Angin')
    ax.set_zlabel('Amplitudo Gelombang')
    ax.set_zlim(-7, 7)

    # Menambahkan anotasi untuk amplitudo dan panjang gelombang
    amplitudo = np.max(height_map)
    panjang_gelombang = np.pi * 2 / np.max(wind_angles)
    ax.text2D(0.05, 0.95, f'Amplitudo: {amplitudo:.2f}', transform=ax.transAxes, color='blue', fontsize=10)
    ax.text2D(0.05, 0.90, f'Panjang Gelombang: {panjang_gelombang:.2f}', transform=ax.transAxes, color='red', fontsize=10)

    return surface,

# Membuat animasi
ani = FuncAnimation(fig, update, frames=range(150), interval=50)
# ani.save('water.mp4', writer='ffmpeg', fps=30, dpi=300)
plt.show()
