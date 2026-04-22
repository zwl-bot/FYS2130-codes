import numpy as np
import matplotlib.pyplot as plt

#from matplotlib.pyplot import cm

#cmap = plt.colormaps["viridis"]


########## Initialization ##########
def EM_wave_arrays(N, y_range, E0, k, omega, t):
    x_arr = np.zeros(N+1)
    y_arr = np.linspace(y_range[0], y_range[1], N+1)
    z_arr = np.zeros(N+1)
    
    
    x_arr = E0*np.cos(k*y_arr-omega*t)
    z_arr = E0*np.sin(k*y_arr-omega*t)  
    
    return x_arr, y_arr, z_arr



########## Plot electric field trajectory ##########
def EM_wave_plot_3d(x_arr, y_arr, z_arr):
    fig = plt.figure(); ax = plt.axes(projection='3d')
    
    
    ax.plot3D(x_arr, y_arr, z_arr, 'blue', label="E(y,t)")
    ax.set_title('3D Electric field')
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.plot3D(np.zeros(len(y_arr)), y_arr, np.zeros(len(y_arr)))
    plt.legend(); plt.show()
    return



########## To verify circular polarization ##########
def EM_wave_plot_2d(x_arr, z_arr, E0):
    plt.plot(0, 0, color="red", marker="o")
    plt.plot([0, E0], [0, 0], color="green", linestyle="dashed", label="E0")
    plt.plot(x_arr, z_arr, label="E(y,t)")
    plt.axis("square")
    plt.title("2D Electric field")
    plt.xlabel("x"); plt.ylabel("z")
    plt.legend(); plt.show()
    return


#---------------------------------------------------------#
                ########## ACTION ##########
#---------------------------------------------------------#
if __name__=="__main__":
    print("HIS 5 - Problem 2")
    
    ##### Initial conditions #####
    y_range = [-10,0]
    E0 = 2
    k = 1
    omega = 1
    t = 1
    N = 100
    
    
    
    ##### Make computations #####
    x_arr, y_arr, z_arr = EM_wave_arrays(N, y_range, E0, k, omega, t)
    
    
    
    ##### Plot time and frequency domains #####
    #EM_wave_plot_3d(x_arr, y_arr, z_arr)
    
    EM_wave_plot_2d(x_arr, z_arr, E0)
    
    
