import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Function to initialize the plot
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    return line,

# Function to update the plot for each frame
def update(frame):
    y = np.sin(2 * np.pi * (x - c * frame / frames))  # Equation for the wave
    line.set_ydata(y)
    return line,

# Parameters
c = 2  # Wave speed
length = 10.0  # Length of the x-axis
frames = 100  # Number of frames in the animation

# Create x-axis values
x = np.linspace(0, length, 1000)

# Set up the plot
fig, ax = plt.subplots()
line, = ax.plot(x, np.sin(2 * np.pi * x), lw=2)  # Initial wave

# Set plot limits
ax.set_ylim(-5, 5)
ax.set_xlim(0, length)

# Create the animation
animation = FuncAnimation(fig, update, frames=np.arange(0, frames), init_func=init, blit=True)

# Display the animation
plt.show()
