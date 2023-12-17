import numpy as np
import matplotlib.pyplot as plt

# Plate dimensions and initial conditions
width = 50  # Number of grid points along the x-axis
height = 50  # Number of grid points along the y-axis
zero = np.zeros((1, 1))

# Initialize the temperature array
T = np.zeros((height, width))

# Set boundary conditions
T[:, 0] = np.linspace(0,120, len(T[0, :]))  # Left side
T[:, -1] = 50  # Right side
T[0, :] = np.linspace(0,0, len(T[0, :]))  # Top side
T[-1, :] = 100   # Bottom side

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
fig, ax = plt.subplots()

im1 = ax.imshow(T, cmap='hot',extent=[0, 1, 0, 1])
# im1 = ax.imshow(image_data, cmap='viridis')
im2 = ax.imshow(zero, cmap='plasma', extent=[0.5, 1, 0.5, 1])
# im3 = ax.imshow(image_data, cmap='hot', extent=[0, 10, 0, 10])

# Set axis limits
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

# Set axis labels
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')

# Add colorbar
cbar = fig.colorbar(im1, ax=ax, orientation='vertical')


# Show the plot
plt.show()
