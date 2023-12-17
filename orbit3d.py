import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as c
import time
import sys
import matplotlib.animation as animation


class Body():
    def __init__(self, mass, x_vec, v_vec, name=None, has_units=True):
        '''
        Membuat instance dari kelas Body, yang digunakan dalam Simulasi.

        :param: mass | massa benda. Jika has_units=True, berupa Astropy Quantity, jika tidak, berupa float.
        :param: x_vec | sebuah vektor len(3) yang berisi posisi awal x, y, z dari benda.
                        Array dapat berupa unitless jika has_units=False, atau berupa np.array([0,0,0])*u.km jika has_units=True.
        :param: v_vec | vektor len(3) yang berisi kecepatan awal v_x, v_y, v_z dari benda.
        :param: name | string yang berisi nama, digunakan untuk plotting nanti.
        :param: has_units | mendefinisikan bagaimana kode memperlakukan masalah, apakah ber-unit atau tidak.
        '''
        self.name = name
        self.has_units = has_units
        if self.has_units:
            self.mass = mass.cgs
            self.x_vec = x_vec.cgs.value
            self.v_vec = v_vec.cgs.value
        else:
            self.mass = mass
            self.x_vec = x_vec
            self.v_vec = v_vec

    def return_vec(self):
        '''
        Menggabungkan vektor x dan v menjadi 1 vektor 'y' yang akan digunakan dalam Runge-Kutta.
        '''
        return np.concatenate((self.x_vec, self.v_vec))

    def return_mass(self):
        '''
        Handler untuk menghilangkan satuan massa jika ada (setelah dikonversi ke cgs) atau mengembalikan nilai float.
        '''
        if self.has_units:
            return self.mass.cgs.value
        else:
            return self.mass

    def return_name(self):
        return self.name


class Simulation():
    def __init__(self, bodies, has_units=True):
        '''
        Inisialisasi instance dari objek Simulasi.
        -------------------------------------------
        Params:
            bodies (list): sebuah list dari objek Body()
            has_units (bool): menentukan apakah bodies yang dimasukkan memiliki satuan atau tidak.
        '''
        self.has_units = has_units
        self.bodies = bodies
        self.N_bodies = len(self.bodies)
        self.nDim = 6.0
        self.quant_vec = np.concatenate(np.array([i.return_vec() for i in self.bodies]))
        self.mass_vec = np.array([i.return_mass() for i in self.bodies])
        self.name_vec = [i.return_name() for i in self.bodies]

    def set_diff_eq(self, calc_diff_eqs, **kwargs):
        '''
        Mengatur fungsi yang menghitung persamaan diferensial yang akan diselesaikan.
        --------------------------------
        Params:
            calc_diff_eqs: fungsi yang menghitung persamaan diferensial.
            kwargs: keyword arguments yang akan digunakan dalam fungsi.
        '''
        self.diff_eq_kwargs = kwargs
        self.calc_diff_eqs = calc_diff_eqs

    def rk4(self, t, dt):
        '''
        Integrator Runge-Kutta orde 4.
        -------------------------------
        Params:
            t: waktu saat ini
            dt: timestep yang akan digunakan
        Returns:
            y_new: vektor y baru yang telah diintegrasi.
        '''
        k1 = dt * self.calc_diff_eqs(t, self.quant_vec, self.mass_vec, **self.diff_eq_kwargs)
        k2 = dt * self.calc_diff_eqs(t + 0.5 * dt, self.quant_vec + 0.5 * k1, self.mass_vec, **self.diff_eq_kwargs)
        k3 = dt * self.calc_diff_eqs(t + 0.5 * dt, self.quant_vec + 0.5 * k2, self.mass_vec, **self.diff_eq_kwargs)
        k4 = dt * self.calc_diff_eqs(t + dt, self.quant_vec + k2, self.mass_vec, **self.diff_eq_kwargs)

        y_new = self.quant_vec + ((k1 + 2 * k2 + 2 * k3 + k4) / 6.0)

        return y_new
    

    def run(self, T, dt, t0=0):
        '''
        Metode yang menjalankan simulasi pada serangkaian benda.
        ---------------------
        Params:
            T: waktu total (dalam satuan simulasi) untuk menjalankan simulasi. Dapat memiliki satuan atau tidak, atur has_units dengan benar.
            dt: timestep (dalam satuan simulasi) untuk memajukan simulasi. Sama seperti di atas.
            t0 (opsional): atur waktu mulai simulasi yang tidak nol.

        Returns:
            None, tetapi menyimpan atribut history yang diakses melalui
            'simulation.history' yang berisi semua vektor y untuk simulasi.
            Ini berbentuk (Nstep, Nbodies * 6), sehingga posisi x dan y partikel 1 adalah
            simulation.history[:,0], simulation.history[:,1], sedangkan yang sama untuk partikel 2 adalah
            simulation.history[:,6], simulation.history[:,7]. Kecepatan juga dapat diekstraksi.
        '''
        if not hasattr(self, 'calc_diff_eqs'):
                raise AttributeError('Anda harus mengatur pemecah persamaan diferensial terlebih dahulu.')
        if self.has_units:
            try:
                _ = t0.unit
            except:
                t0 = (t0 * T.unit).cgs.value
            T = T.cgs.value
            dt = dt.cgs.value

        self.history = [self.quant_vec]
        clock_time = t0
        nsteps = int((T - t0) / dt)
        start_time = time.time()
        for step in range(nsteps):
            sys.stdout.flush()
            sys.stdout.write('Mengintegrasi: langkah = {} / {} | waktu simulasi = {}'.format(step, nsteps, round(clock_time, 3)) + '\r')
            y_new = self.rk4(0, dt)
            self.history.append(y_new)
            self.quant_vec = y_new
            clock_time += dt
        runtime = time.time() - start_time
        print('\n')
        print('Simulasi selesai dalam {} detik'.format(runtime))
        self.history = np.array(self.history)

def nbody_solve(t, y, masses):
    N_bodies = int(len(y) / 6)
    solved_vector = np.zeros(y.size)
    for i in range(N_bodies):
        ioffset = i * 6
        for j in range(N_bodies):
            joffset = j * 6
            solved_vector[ioffset] = y[ioffset + 3]
            solved_vector[ioffset + 1] = y[ioffset + 4]
            solved_vector[ioffset + 2] = y[ioffset + 5]
            if i != j:
                dx = y[ioffset] - y[joffset]
                dy = y[ioffset + 1] - y[joffset + 1]
                dz = y[ioffset + 2] - y[joffset + 2]
                r = (dx ** 2 + dy ** 2 + dz ** 2) ** 0.5
                ax = (-c.G.cgs * masses[j] / r ** 3) * dx
                ay = (-c.G.cgs * masses[j] / r ** 3) * dy
                az = (-c.G.cgs * masses[j] / r ** 3) * dz
                ax = ax.value
                ay = ay.value
                az = az.value
                solved_vector[ioffset + 3] += ax
                solved_vector[ioffset + 4] += ay
                solved_vector[ioffset + 5] += az
    return solved_vector

M_sun = (1.989e30 * u.kg).cgs  # Massa Matahari dalam satuan CGS

Sun = Body(mass=M_sun,
           x_vec=np.array([0, 0, 0]) * u.km,
           v_vec=np.array([0, 0, 0]) * u.km / u.s,
           name='Matahari')

v_earth = np.array([0, 29.78, 0]) * u.km / u.s  # Kecepatan Bumi relatif terhadap Matahari

Earth = Body(mass=c.M_earth.cgs,
             x_vec=np.array([1.496e8, 0, 0]) * u.km,  # Jarak Bumi dari Matahari dalam kilometer
             v_vec=v_earth,
             name='Bumi')

bodies = [Sun, Earth]

simulation = Simulation(bodies)
simulation.set_diff_eq(nbody_solve)

simulation.run(366*u.day, 5*u.hr)

fig = plt.figure(figsize=(8, 8), tight_layout=True)
ax = fig.add_subplot(111, projection='3d')

# Plotting the initial positions of Earth and Sun
sun_dot, = ax.plot(simulation.history[0, 0], simulation.history[0, 1], 'o', color='gold', markersize=12, label='Matahari', zorder=10, alpha=0.8)
earth_dot, = ax.plot(simulation.history[0, 6], simulation.history[0, 7], 'o', color='teal', markersize = 7, label='Bumi', zorder=10, alpha=0.8)

# Plotting the orbit lines of Earth and Sun
earth_orbit, = ax.plot(simulation.history[:, 6], simulation.history[:, 7],  label='Orbit Bumi', alpha=0.5)
sun_orbit, = ax.plot(simulation.history[:, 0], simulation.history[:, 1], alpha=0.5)

ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_zlabel('z (km)')
ax.legend()

def update(frame):
    # Update the positions of Earth and Moon dots
    earth_dot.set_data(simulation.history[frame, 6], simulation.history[frame, 7])
    earth_dot.set_3d_properties(0)
    
    return earth_dot,

ani = animation.FuncAnimation(fig, update, frames=len(simulation.history), interval=50, blit=True)
# ani.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'], dpi=300)
plt.show()