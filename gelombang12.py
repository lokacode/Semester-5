import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parameters
speed = 0.5  # Wave packet speed
x_max = 20   # Maximum x-axis value
num_points = 500  # Number of points along the x-axis
# Randomize barrier position and height
# num_barriers = np.random.randint(20, 24)
num_barriers = 12
barrier_positions = np.sort(np.linspace(7.5, 12.5, num_barriers))
barrier_heights = np.sort(np.random.uniform(0.05, 0.1, num_barriers))
barrier_heights = barrier_heights[::-1]

print('Number of barriers:', num_barriers)
print('Barrier positions:', barrier_positions)
print('Barrier heights:', barrier_heights)

barrier_data = pd.DataFrame({
    'Barrier_Positions': barrier_positions,
    'Barrier_Heights': barrier_heights
})

# Save the DataFrame to an Excel file
barrier_data.to_excel('barrier_data15.xlsx', index=False)

# Create x-axis values
x_values = np.linspace(0, x_max, num_points)

# Create a function for the continuous sine wave with a potential barrier
def sine_wave_with_barrier(x, t):
    wave = np.sin(10*np.pi * x / x_max - speed * t)
    
    # Calculate the potential barrier
    barrier = 0
    for i in range(num_barriers):
        barrier += np.where(x > barrier_positions[i], barrier_heights[i], 0)

    # Add the potential barrier
    wave *= np.exp(-barrier)

    return wave

# Create an animation
fig, ax = plt.subplots()
line, = ax.plot(x_values, sine_wave_with_barrier(x_values, 0), color='blue', label='Gelombang air')
barrier_lines = [ax.axvline(pos, color='forestgreen', lw=1) for pos in barrier_positions]
ax.axvline(barrier_positions[0], color='forestgreen', lw=1, label='Pohon mangrove')
ax.axvspan(7.5, 12.5, color='forestgreen', alpha=0.2, lw=0, label='Area mangrove')
ax.legend(loc='upper right')

def update(frame):
    line.set_ydata(sine_wave_with_barrier(x_values, frame))
    return line,

ani = FuncAnimation(fig, update, frames=np.arange(0, x_max/speed, 0.05), interval=5, blit=True)

# Set up plot
ax.set_xlim(0, x_max)
ax.set_ylim(-5, 5)
ax.set_yticks(np.arange(-5, 6))
ax.set_xlabel('Jarak (m)')
ax.set_title(f'Pengaruh pohon mangrove terhadap gelombang air (Jumlah pohon: {num_barriers})')
ax.hlines(0, 0, x_max, color='black', lw=1)
ax.set_ylabel('Amplitudo gelombang (m)')
plt.tight_layout()
ani.save('gelombang15.mp4', dpi=300, fps=30, writer='ffmpeg')
plt.show()