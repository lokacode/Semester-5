import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
import astropy.constants as c 

class Body():
    def __init__(self, name, x_vec, v_vec, mass):
        self.name = name
        self.x_vec = x_vec
        self.v_vec = v_vec
        self.a_vec = np.array([0,0,0])*u.km/u.s/u.s
        self.mass = mass
    
    def return_vec(self):
        return self.v_vec
    

Earth = Body(name='Earth',
             x_vec = np.array([0,1,0])*u.AU,
             v_vec = np.array([0,30,0])*u.km/u.s,
             mass = 1.0*c.M_earth) 
Sun = Body(name='Sun',
           x_vec = np.array([0,0,0])*u.AU,
           v_vec = np.array([0,0,0])*u.km/u.s,
           mass = 1*u.Msun)

bodies = [Earth,Sun]

def set_diff_eq(self,calc_diff_eqs,**kwargs):
    '''
    Method which assigns an external solver function as the diff-eq solver for RK4. 
    For N-body or gravitational setups, this is the function which calculates accelerations.
    ---------------------------------
    Params:
        calc_diff_eqs: A function which returns a [y] vector for RK4
        **kwargs: Any additional inputs/hyperparameters the external function requires
    '''
    self.calc_diff_eqs = calc_diff_eqs
    self.kwargs = kwargs

def rk4(t,dt,y,evaluate):
    '''
    Given a vector y with [x,xdot], calculate the new y for step dt,
    using rk4 method
    '''
    k1 = dt * evaluate(t, y) 
    k2 = dt * evaluate(t + 0.5*dt, y + 0.5*k1)
    k3 = dt * evaluate(t + 0.5*dt, y + 0.5*k2)
    k4 = dt * evaluate(t + dt, y + k3)
    
    y_new = y + (1/6.)*(k1+ 2*k2 + 2*k3 + k4)
    return y_new

def evaluate_SHO(t,y,k=1):
    '''
    evaluate the SHO at time t and y=y. 
    Note: this diff eq has no time dependence
    '''
    v = y[1] 
    a = -k**2 * y[0]
    return np.array([v,a])

# Running a small integration

y_0 = np.array([-5,0]) #initialize oscillator at x = -5, with 0 velocity. 
history = [y_0]
ts = [0]
dt = 0.01
T = 10
nsteps = int(T/dt)
for i in range(nsteps):
    y_new = rk4(T,dt,history[-1],evaluate_SHO)
    history.append(y_new)
    t = ts[-1] + dt
    ts.append(t)
history = np.array(history)
ts = np.array(ts)

analytical_solution = -5*np.cos(ts)
analytical_velocity = 5*np.sin(ts)
fig, ax = plt.subplots(2,1,figsize=(12,6),sharex=True)
ax[0].plot(ts,history[:,0],color='C0',lw=6,ls='--',label='Position (rk4)',alpha=0.5)
ax[0].plot(ts,analytical_solution,color='r',label='Analytical Solution')
ax[1].plot(ts,history[:,1],color='C1',lw=6,alpha=0.5,ls='--',label='Velocity (rk4)')
ax[1].plot(ts,analytical_velocity,'C2',label='Analytical Solution')
ax[0].legend(loc="upper center")
ax[1].legend(loc="upper center")
ax[-1].set_xlabel('time')
# plt.show()

class Simulation():
    def __init__(self,bodies):
        self.bodies = bodies
        self.N_bodies = len(bodies)
        self.Ndim = 6
        self.quant_vec = np.zeros(self.N_bodies*self.Ndim)

    def set_diff_eq(self,calc_diff_eqs,**kwargs):
        '''
        Method which assigns an external solver function as the diff-eq solver for RK4. 
        For N-body or gravitational setups, this is the function which calculates accelerations.
        ---------------------------------
        Params:
            calc_diff_eqs: A function which returns a [y] vector for RK4
            **kwargs: Any additional inputs/hyperparameters the external function requires
        '''
        self.calc_diff_eqs = calc_diff_eqs
        self.kwargs = kwargs

    def rk4(self,t,dt):
        '''
        RK4 integrator. Calculates the K values and returns a new y vector
        --------------------------------
        Params:
            t: a time. Only used if the diff eq depends on time (gravity doesn't).
            dt: timestep. Non adaptive in this case
        '''
        k1 = dt * self.calc_diff_eqs(t, self.quant_vec,self.mass_vec,**self.diff_eq_kwargs) 
        k2 = dt * self.calc_diff_eqs(t + 0.5*dt, self.quant_vec + 0.5*k1, self.mass_vec, **self.diff_eq_kwargs)
        k3 = dt * self.calc_diff_eqs(t + 0.5*dt, self.quant_vec + 0.5*k2, self.mass_vec, **self.diff_eq_kwargs)
        k4 = dt * self.calc_diff_eqs(t + dt, self.quant_vec + k3, self.mass_vec, **self.diff_eq_kwargs)
        y_new = self.quant_vec + ((k1 + 2*k2 + 2*k3 + k4) / 6.0)

        return y_new
    
    def run(self,T,dt,t0=0):
        '''
        Method which runs the simulation on a given set of bodies.
        ---------------------
        Params: 
            T: total time (in simulation units) to run the simulation. Can have units or not, just set has_units appropriately.
            dt: timestep (in simulation units) to advance the simulation. Same as above
            t0 (optional): set a non-zero start time to the simulation.

        Returns: 
            None, but leaves an attribute history accessed via 
            'simulation.history' which contains all y vectors for the simulation. 
            These are of shape (Nstep,Nbodies * 6), so the x and y positions of particle 1 are
            simulation.history[:,0], simulation.history[:,1], while the same for particle 2 are
            simulation.history[:,6], simulation.history[:,7]. Velocities are also extractable.
        '''
        pass

def nbody_solve(t,y,masses):
    N_bodies = int(len(y) / 6)
    solved_vector = np.zeros(y.size)
    for i in range(N_bodies):
        ioffset = i * 6 
        for j in range(N_bodies):
            joffset = j*6
            solved_vector[ioffset] = y[ioffset+3]
            solved_vector[ioffset+1] = y[ioffset+4]
            solved_vector[ioffset+2] = y[ioffset+5]
            if i != j:
                dx = y[ioffset] - y[joffset]
                dy = y[ioffset+1] - y[joffset+1]
                dz = y[ioffset+2] - y[joffset+2] 
                r = (dx**2+dy**2+dz**2)**0.5
                ax = (-c.G.cgs * masses[j] / r**3) * dx
                ay = (-c.G.cgs * masses[j] / r**3) * dy
                az = (-c.G.cgs * masses[j] / r**3) * dz
                ax = ax.value
                ay = ay.value
                az = az.value
                solved_vector[ioffset+3] += ax
                solved_vector[ioffset+4] += ay
                solved_vector[ioffset+5] += az            
    return solved_vector 

Earth = Body(mass=c.M_earth.cgs,
             x_vec=np.array([0,0,0])*u.km,
             v_vec=np.array([0,0,0])*u.km/u.s,
             name='Earth')


M_moon = (7.347e22*u.kg).cgs
Moon = Body(mass=M_moon,
           x_vec = np.array([3.84e5,0,0])*u.km,
           v_vec = np.array([0,1.022,0])*u.km/u.s,
           name='Moon')


bodies = [Earth,Moon]

simulation = Simulation(bodies)
simulation.set_diff_eq(nbody_solve)

simulation.run(72*u.day,1*u.hr)

# Extracting the positions of Earth and Moon from the history
earth_positions = history[:, 0:3]
moon_positions = history[:, 6:9]

# Plotting the positions of Earth and Moon
fig, ax = plt.subplots(figsize=(6, 6))

# Plotting the initial positions of Earth and Sun
sun_dot, = ax.plot(simulation.history[0, 0], simulation.history[0, 1], 'yo', markersize = 12, label='Matahari', zorder=10)
earth_dot, = ax.plot(simulation.history[0, 6], simulation.history[0, 7], 'bo', markersize = 7, label='Bumi', zorder=10)

# Plotting the orbit lines of Earth and Sun
earth_orbit, = ax.plot(simulation.history[:, 6], simulation.history[:, 7],  label='Orbit Bumi')
Sun_orbit, = ax.plot(simulation.history[:, 0], simulation.history[:, 1])

ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_xlim(-1.7e13, 1.7e13)
ax.set_ylim(-1.7e13, 1.7e13)
ax.set_title('Simulasi Orbit Bumi dan Bulan dengan RK4')
ax.legend()

def update(frame):
    # Update the positions of Earth dots
    earth_dot.set_xdata(simulation.history[frame, 6])
    earth_dot.set_ydata(simulation.history[frame, 7])
    
    return earth_dot

ani = animation.FuncAnimation(fig, update, frames=len(simulation.history), interval=10, blit=True)
ani.save('orbit.mp4', fps=30, extra_args=['-vcodec', 'libx264'], dpi=300)
# plt.show()