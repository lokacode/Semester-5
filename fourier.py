import numpy as np
import matplotlib.pyplot as plt

omega = 2 * np.pi  # Frekuensi sinyal

# f(x) = 4/pi * (sin(omega * x) + 1/3 * sin(3 * omega * x) + 1/5 * sin(5 * omega * x) + ... + 1/n * sin(n * omega * x))
def fourier_series1(x, n):
    result = 0
    for i in range(1, n + 1):
        result += (1 / (2 * i - 1)) * np.sin((2 * i - 1) * omega * x)
    result *= 4 / np.pi
    return result


# f(x) = 3/2 + 6/pi * (sin(pi * x / 5) + 1/3 * sin(3 * pi * x / 5) + 1/5 * sin(5 * pi * x / 5) + ... + 1/n * sin(n * pi * x / 5))
def fourier_series2(x, n):
    result = 3/2
    for i in range(1, n + 1):
        result += (6 / np.pi) * (1 / (2 * i - 1)) * np.sin((2 * i - 1) * np.pi * x / 5)
    return result

x1 = np.linspace(0, np.pi, 1000)
x2 = np.linspace(-10, 10, 1000)

n = 10
y1 = fourier_series1(x1, n)
y2 = fourier_series2(x2, n)

# export data to xlsx
# import pandas as pd
# df = pd.DataFrame({'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2})
# df.to_excel('fourier.xlsx', sheet_name='Sheet1', index=False)

ax = plt.subplot(2, 1, 1)
ax.plot(x1, y1, 'r')
ax.set_title(f'Tugas 1')
ax.set_xlabel('$x$')
ax.set_ylabel(f'$y$')
ax.grid()

ax = plt.subplot(2, 1, 2)
ax.plot(x2, y2,'b')
ax.set_title(f'Tugas 2')
ax.set_xlabel('$x$')
ax.set_ylabel(f'$y$')
ax.grid()

plt.tight_layout()
# plt.savefig('fourier.png', dpi=300)
plt.show()