import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parameters
size = 100  # Size of the grid
dx = 1.0    # Grid spacing
dt = 0.1    # Time step
c = 1.0     # Wave speed

# Create grid
x = np.arange(0, size) * dx
y = np.arange(0, size) * dx
X, Y = np.meshgrid(x, y)

# Initial wave function (you can modify this to create different initial conditions)
initial_wave = np.exp(-((X - size/2)**2 + (Y - size/2)**2) / 10)

# Initialize wave arrays
current_wave = initial_wave.copy()
next_wave = np.zeros_like(current_wave)

# Function to update the wave over time
def update_wave(frame):
    global current_wave, next_wave

    # Compute the wave using the wave equation
    next_wave[1:-1, 1:-1] = 2 * current_wave[1:-1, 1:-1] - next_wave[1:-1, 1:-1] + \
                            (c**2 * dt**2 / dx**2) * (
                                current_wave[2:, 1:-1] +
                                current_wave[:-2, 1:-1] +
                                current_wave[1:-1, 2:] +
                                current_wave[1:-1, :-2] -
                                4 * current_wave[1:-1, 1:-1]
                            )

    # Update current wave for the next iteration
    current_wave, next_wave = next_wave, current_wave

    # Update the plot with the new wave
    ax.clear()
    ax.imshow(current_wave, cmap='viridis', origin='lower')

# Set up the plot
fig, ax = plt.subplots()
ax.imshow(initial_wave, cmap='viridis', origin='lower')
ax.set_title('Wave Equation')

# Create animation
animation = FuncAnimation(fig, update_wave, frames=200, interval=50)

# save animation
# animation.save('wave.mp4', writer='ffmpeg', fps=30, dpi=300)

# Show the animation
plt.show()
