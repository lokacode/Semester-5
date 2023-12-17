import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
air_density = 1.225  # kg/m^3
water_density = 1000  # kg/m^3
iron_density = 7874  # kg/m^3

air_speed = 343  # m/s (speed of sound in air at room temperature)
water_speed = 1500  # m/s (speed of sound in water)
iron_speed = 5000  # m/s (speed of sound in iron)

# Time parameters
total_time = 2.0  # seconds
dt = 0.01  # time step

# Spatial parameters
distance_air = 50  # meters
distance_water = 20  # meters
distance_iron = 10  # meters

# Create time array
time_points = np.arange(0, total_time, dt)

# Create distance arrays for each medium
distance_air_array = np.linspace(0, distance_air, len(time_points))
distance_water_array = np.linspace(distance_air, distance_air + distance_water, len(time_points))
distance_iron_array = np.linspace(distance_air + distance_water, distance_air + distance_water + distance_iron, len(time_points))

# Create speed arrays for each medium
speed_air = np.full_like(distance_air_array, air_speed)
speed_water = np.full_like(distance_water_array, water_speed)
speed_iron = np.full_like(distance_iron_array, iron_speed)

# Concatenate arrays
distance_array = np.concatenate((distance_air_array, distance_water_array, distance_iron_array))
speed_array = np.concatenate((speed_air, speed_water, speed_iron))

# Initialize plot
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, distance_air + distance_water + distance_iron)
ax.set_ylim(0, max(speed_air.max(), speed_water.max(), speed_iron.max()) + 1000)
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Speed of Sound (m/s)')

# Animation function
def update(frame):
    line.set_data(distance_array[:frame], speed_array[:frame])
    return line,

# Create animation
animation = FuncAnimation(fig, update, frames=len(distance_array), interval=dt * 1000, blit=True)

plt.show()
