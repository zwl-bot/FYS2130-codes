import numpy as np
import matplotlib.pyplot as plt


#---------------------------------------------------------#
            ########## d) (iii) ##########
#---------------------------------------------------------#


##### wave function u(x, t; k, w) #####
def wave_func(x, t, k, w): return np.cos(k*x - w*t)



##### Create Wave function 2D array #####
def wave_func_2D(x_interval, T, dx, dt, k, w):
    x_arr = np.arange(x_interval[0], x_interval[1] + dx, dx)
    t_arr = np.arange(0, T + dt, dt)
    u = np.ones((len(t_arr), len(x_arr)))
    
    for t in range(len(t_arr)):
        for x in range(len(x_arr)):
            u[t, x] = wave_func(x_arr[x], t_arr[t], k, w)
    return u, x_arr, t_arr



##### Plot wave equation at different times #####
def wave_func_2D_plot_snapshots(u_arr, x_arr, t_arr):
    snap_step = int(len(t_arr)/4)
    t_snaps = [0, snap_step, snap_step*2, snap_step*3, snap_step*4]
    for snap in t_snaps:
        plt.plot(x_arr, u_arr[snap], label=f"t={t_arr[snap]:.1f}s")
    plt.legend(); plt.title("Wave equation different time snapshots")
    plt.xlabel("x"); plt.ylabel("u(x,t)"); plt.show()





#---------------------------------------------------------#
                ########## ACTION ##########
#---------------------------------------------------------#
if __name__ == "__main__":
    print("Problem 2 d) (iii)")
    
    ##### Initial conditions #####
    interval = [0, 1]; dx = 0.01        # interval x-axis; x-step dx
    T = 5; dt = 0.1                     # Time span T (s); Timestep dt
    k = 10                              # Wavenumber k
    omega = 20                          # Frequency w
    
    
    ##### Make computations #####
    u, xs, ts = wave_func_2D(interval, T, dx, dt, k, omega)
    
    
    ##### Plot snapshot at different times #####
    wave_func_2D_plot_snapshots(u, xs, ts)
    