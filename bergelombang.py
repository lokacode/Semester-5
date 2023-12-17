import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
air_speed = 343  # speed of sound in air (m/s)
water_speed = 1481  # speed of sound in water (m/s)
iron_speed = 5130  # speed of sound in iron (m/s)
f = 1/1000
def wave(x):
    return np.sin(f * np.pi * x)

# x = np.linspace(0, 10, 1000)
# y = wave(x)
# plt.plot(x, y)
# Parameters
duration_air = 2  # duration of air simulation in seconds
duration_water = 2  # duration of water simulation in seconds
duration_iron = 2  # duration of iron simulation in seconds
fps = 30  # frames per second

# Calculate distances
distance_air = air_speed * duration_air
distance_water = water_speed * duration_water
distance_iron = iron_speed * duration_iron

# Create time arrays
time_air = np.linspace(0, duration_air, int(fps * duration_air))
time_water = np.linspace(duration_air, duration_air + duration_water, int(fps * duration_water))
time_iron = np.linspace(duration_air + duration_water, duration_air + duration_water + duration_iron, int(fps * duration_iron))

a = np.linspace(0, distance_air, len(time_air))
b = np.linspace(distance_air, distance_air + distance_water, len(time_water))
c = np.linspace(distance_air + distance_water, distance_air + distance_water + distance_iron, len(time_iron))
print(len(a), len(b), len(c))
# Create distance arrays
distance = np.concatenate([a,
                           b,
                           c])

distance2 = np.linspace(0, max(distance), len(distance))
# Create speed arrays
speed = np.concatenate([np.full_like(time_air, air_speed),
                        np.full_like(time_water, water_speed),
                        np.full_like(time_iron, iron_speed)])

speed = np.linspace(0, 0, len(distance))

# Create the plot
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, distance_iron+distance_air+distance_water)
ax.set_ylim(-iron_speed-500, iron_speed + 500)
ax.set_xlabel('Distance (m)')
ax.axvline(distance_air, color="teal")
ax.axvline(distance_air+distance_water, color="blue")
ax.set_ylabel('Speed of Sound (m/s)')

# Update function for animation
def update(frame):
    line.set_data(distance[:frame], (speed[:frame] + wave(distance2[:frame]) * 1000))
    return line,

# Create animation
animation = FuncAnimation(fig, update, frames=len(distance), interval=1000 / fps)

plt.show()
