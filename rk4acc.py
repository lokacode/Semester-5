import numpy as np

def rk4(t,dt,y,evaluate):
    '''
    Diberikan vektor y dengan [x,xdot], hitung y baru untuk langkah dt,
    menggunakan metode rk4
    '''
    k1 = dt * evaluate(t, y) 
    k2 = dt * evaluate(t + 0.5*dt, y + 0.5*k1)
    k3 = dt * evaluate(t + 0.5*dt, y + 0.5*k2)
    k4 = dt * evaluate(t + dt, y + k3)
    
    y_new = y + (1/6.)*(k1+ 2*k2 + 2*k3 + k4)
    return y_new

def evaluate_SHO(t,y,k=1):
    '''
    evaluasi SHO pada waktu t dan y=y. 
    Catatan: persamaan diferensial ini tidak memiliki ketergantungan waktu
    '''
    v = y[1] 
    a = -k**2 * y[0]
    return np.array([v,a])

# Menjalankan integrasi kecil

y_0 = np.array([-5,0]) #inisialisasi osilator pada x = -5, dengan kecepatan 0. 
history = [y_0]
ts = [0]
dt = 0.01
T = 17.5
nsteps = int(T/dt)
for i in range(nsteps):
    y_new = rk4(T,dt,history[-1],evaluate_SHO)
    history.append(y_new)
    t = ts[-1] + dt
    ts.append(t)
history = np.array(history)
ts = np.array(ts)


import matplotlib.pyplot as plt

solusi_analitis = -5*np.cos(ts)
kecepatan_analitis = 5*np.sin(ts)

position_error = np.abs(history[:, 0] - solusi_analitis)
velocity_error = np.abs(history[:, 1] - kecepatan_analitis)

position_rmse = np.sqrt(np.mean(position_error**2))
velocity_rmse = np.sqrt(np.mean(velocity_error**2))

print(f"Position RMSE: {position_rmse}")
print(f"Velocity RMSE: {velocity_rmse}")

fig, ax = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
ax[0].plot(ts, position_error, color='C0', lw=2, label='Position Error')
ax[1].plot(ts, velocity_error, color='C1', lw=2, label='Velocity Error')
ax[0].set_ylabel('Position Error')
ax[1].set_ylabel('Velocity Error')
ax[-1].set_xlabel('Time')
ax[0].legend(loc="upper left")
ax[1].legend(loc="upper left")
# plt.savefig('rk4SHO_error.png', dpi=300)
plt.show()

# fig, ax = plt.subplots(2,1,figsize=(12,6),sharex=True)
# ax[0].plot(ts,history[:,0],color='C0',lw=6,ls='--',label='Posisi (rk4)',alpha=0.5)
# ax[0].plot(ts,solusi_analitis,color='r',label='Solusi Analitis')
# ax[1].plot(ts,history[:,1],color='C1',lw=6,alpha=0.5,ls='--',label='Kecepatan (rk4)')
# ax[1].plot(ts,kecepatan_analitis,'C2',label='Solusi Analitis')
# ax[0].legend(loc="upper center")
# ax[1].legend(loc="upper center")
# ax[-1].set_xlabel('waktu')
# ax[0].set_title('Perbandingan numerik dan analitis untuk mengukur ketepatan simulasi')
# plt.savefig('rk4SHO.png',dpi=300)
# plt.show()