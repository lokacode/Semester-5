import matplotlib.pyplot as plt

# Definisikan fungsi untuk mengubah posisi titik dan menampilkan teks
def ubah_titik(x1, y1, nama_titik1, x2, y2, nama_titik2, title):
    plt.cla()
    
    # Plot data set 1
    plt.scatter(x1, y1, s=10, label='Debit 500 LPM')
    for i in range(len(x1)):
        plt.text(x1[i], y1[i], f"({nama_titik1[i]})", fontsize=10)

    # Plot data set 2
    plt.scatter(x2, y2, s=10, label='Debit 300 LPM')
    for i in range(len(x2)):
        plt.text(x2[i], y2[i], f"({nama_titik2[i]})", fontsize=10)

    plt.plot(x1, y1)
    plt.plot(x2, y2)

    plt.title(title)
    plt.xlabel("frekuensi (Hz)")
    plt.ylabel("Kedalaman Air (cm)")
    plt.legend()  # Add legend to distinguish between the two data sets
    plt.show()

# Debit 500 LPM jarak 1 cm
x1 = [12.8, 11.8, 11.2, 10.6, 9.8, 9.4]
y1 = [3, 4, 5, 6, 7, 8]
nama_titik1 = ["N2(3,1)", "N2(4,1)", "N2(5,1)", "N2(6,1)", "N2(7,1)", "N2(8,1)"]

# Debir 300 LPM jarak 1 cm
x2 = [10.2, 9.2, 8.8, 8.6, 8.2, 7.6]
y2 = [3, 4, 5, 6, 7, 8]
nama_titik2 = ["N1(3,1)", "N1(4,1)", "N1(5,1)", "N1(6,1)", "N1(7,1)", "N1(8,1)"]

# Tampilkan grafik
ubah_titik(x1, y1, nama_titik1, x2, y2, nama_titik2, "Debit 500 LPM dan 300 LPM Jarak 1 cm")
