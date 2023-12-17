import numpy as np
import matplotlib.pyplot as plt

# Plate dimensions and initial conditions
width = 50  # Number of grid points along the x-axis
height = 50  # Number of grid points along the y-axis

# Initialize the temperature array
T = np.zeros((height, width))

# Set boundary conditions
T[:, 0] = 75  # Left side
T[:, -1] = 50  # Right side
T[0, :] = np.linspace(0,100, len(T[0, :]))  # Bottom side
T[-1, :] = 100   # Top side

# Convergence criteria
tolerance = 1e-4
max_iterations = 20000

# Gauss-Seidel (Liebmann's method) iteration
for iteration in range(max_iterations):
    for i in range(1, height - 1):
        for j in range(1, width - 1):
            T[i, j] = 0.25 * (T[i - 1, j] + T[i + 1, j] + T[i, j - 1] + T[i, j + 1])

    # Check for convergence
    if np.max(np.abs(T - np.roll(T, 1, axis=0))) < tolerance:
        break

# Plot the temperature distribution
plt.imshow(T, cmap='hot', origin='lower', extent=[0, width, 0, height], interpolation='bicubic')
plt.colorbar(label='Temperature (°C)')
plt.title('Temperature Distribution on Heated Plate')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show()
