# Arden Ahmad
# 121110018

import numpy as np
import matplotlib.pyplot as plt

omega = 2 * np.pi  # Frekuensi sinyal

# f(x) = 4/pi * (sin(omega * x) + 1/3 * sin(3 * omega * x) + 1/5 * sin(5 * omega * x) + ... + 1/n * sin(n * omega * x))
def fourier_series(x, n):
    result = 0
    for i in range(1, n + 1):
        result += (1 / (2 * i - 1)) * np.sin((2 * i - 1) * omega * x)
    result *= 4 / np.pi
    return result

# Membuat nilai x
x = np.linspace(0, 2 * np.pi, 1000)

# Tentukan nilai n untuk membuat plot
n_values = [1, 10, 50]
colors = ['r', 'g', 'b']

# Membuat subplot
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Membuat plot untuk setiap nilai n
for i, (n, color) in enumerate(zip(n_values, colors)):
    y = fourier_series(x, n)
    axes[i].plot(x, y, color=color)
    axes[i].set_title(f'Grafik $y$ terhadap $x$ ($n$={n})')
    axes[i].set_xlabel('$x$')
    axes[i].set_ylabel(f'$y$')
    axes[i].grid()

# Menampilkan plot
plt.tight_layout()
# plt.savefig('tugas 1 fourier.png', dpi=300)
plt.show()
