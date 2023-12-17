import numpy as np
import matplotlib.pyplot as plt

# Fungsi yang akan dicari akarnya
def f(x):
    return np.exp(-x) + 3 - x**2

# Membuat rentang nilai x
x = np.linspace(-2, 4, 100)
y = f(x)

# Inisialisasi tebakan awal dan selang
x0 = 0.5
x1 = 3
delta_x = 0.5

# Iterasi untuk membuat tabel
table = []
while abs(x1 - x0) > 1e-6:
    x0 = x1
    x1 = x0 - f(x0) * delta_x / (f(x0 + delta_x) - f(x0))
    table.append([x0, f(x0), x1])

# Menampilkan tabel
print("Tabel Iterasi:")
print("{:<10} {:<15} {:<10}".format("x0", "f(x0)", "x1"))
for row in table:
    print("{:<10} {:<15} {:<10}".format(row[0], row[1], row[2]))

# menyimpan tabel dalam file xlsx
import pandas as pd
df = pd.DataFrame(table, columns=['x0', 'f(x0)', 'x1'])
df.to_excel('tabel.xlsx', index=False)

# Plot fungsi
plt.plot(x, y, label='f(x)')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.grid(True)
plt.legend()
plt.title('Grafik f(x)')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)

# Plot akar yang ditemukan
plt.scatter(x1, f(x1), color='red', label='Akar', marker='x', zorder = 3)
# tulis nilai akar di grafik
plt.text(x1, f(x1), '  {:.6f}'.format(x1), ha='left', va='bottom', color='red', zorder = 3)
plt.legend()

# Menampilkan grafik
# plt.savefig('grafik fungsi.png', dpi=300)
plt.show()
